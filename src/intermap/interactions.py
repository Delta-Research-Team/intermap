# Created by rglez at 12/10/24

import numba.types
import numpy as np
from numba import njit
from numba.typed import Dict as numba_dict
from numba_kdtree import KDTree as nckd

import intermap.cutoffs as cf


@njit
def create_inter_ids():
    # Create a Numba Dict with specific types
    inter_ids = numba_dict.empty(numba.types.unicode_type, numba.types.int64)

    # undirected 2p interactions
    inter_ids['close_contacts'] = 0
    inter_ids['vdw_contacts'] = 1
    inter_ids['hydrophobic'] = 2

    # directed 2p interactions
    inter_ids['anionic'] = 3
    inter_ids['cationic'] = 4
    inter_ids['metal_donor'] = 5
    inter_ids['metal_acceptor'] = 6

    # directed 3p interactions
    inter_ids['hb_acceptor'] = 7
    inter_ids['hb_donor'] = 8
    return inter_ids


# =============================================================================
# Helper functions
# =============================================================================
@njit(parallel=False)
def calc_dist(d, a):
    """
    Computes the Euclidean distance between two atoms in a molecule

    Args:
        d (ndarray): Coordinates of the first atom (n, 3).
        a (ndarray): Coordinates of the second atom (n, 3).

    Returns:
        float: the Euclidean distance between the two atoms
    """
    n = d.shape[0]
    distances = np.empty(n)

    for i in range(n):
        dx = d[i][0] - a[i][0]
        dy = d[i][1] - a[i][1]
        dz = d[i][2] - a[i][2]
        distances[i] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return distances


@njit(parallel=False)
def calc_angle(d, h, a):
    """
    Computes the angles between sets of three atoms.

    Args:
        d (ndarray): Coordinates of the donor atoms (n, 3).
        h (ndarray): Coordinates of the hydrogen atoms (n, 3).
        a (ndarray): Coordinates of the acceptor atoms (n, 3).

    Returns:
        angle_deg: Angles in degrees for each set of atoms (n,).
    """
    n = d.shape[0]  # Number of triplets
    angle_deg = np.empty(n)  # Array to hold the result angles

    for i in range(n):
        # Compute vectors
        dh = d[i] - h[i]
        ah = a[i] - h[i]

        # Compute dot product and norms
        dot_product = dh[0] * ah[0] + dh[1] * ah[1] + dh[2] * ah[2]
        dh_norm = np.sqrt(dh[0] ** 2 + dh[1] ** 2 + dh[2] ** 2)
        ah_norm = np.sqrt(ah[0] ** 2 + ah[1] ** 2 + ah[2] ** 2)

        # Compute angle
        angle_rad = np.arccos(dot_product / (dh_norm * ah_norm))
        angle_deg[i] = np.rad2deg(angle_rad)  # Convert radians to degrees

    return angle_deg


# =============================================================================
# Interaction functions
# =============================================================================


@njit(parallel=False)
def close_contacts(xyz, k, contact_id,
                   s1_indices, s2_indices,
                   dist_cut, vdw_radii=None):
    """
    Find the close contacts between two selections.

    Args:
        xyz (ndarray): Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id (int): Identifier for the contact type.

        s1_indices (ndarray): Indices of the atoms in s1.
        s2_indices (ndarray): Indices of the atoms in s2.

        dist_cut (float): Cutoff distance for the tree query.
        vdw_radii (ndarray): Van der Waals radii of the atoms.

    Returns:
        contact_indices (list): List of labels for each close contact.
    """

    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], dist_cut)

    # Find the close contacts
    n_contacts = sum([len(x) for x in ball_1])
    contact_indices = np.zeros((n_contacts, 4), dtype=np.int32)

    # Fill the labels
    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            contact_indices[counter][0] = X1
            contact_indices[counter][1] = X2
            contact_indices[counter][2] = k
            contact_indices[counter][3] = contact_id

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = contact_indices[:, 0] == contact_indices[:, 1]
    contact_indices = contact_indices[~idems]

    if vdw_radii is None:
        return contact_indices

    # Discard pairs whose distance is greater than the sum of the VdW radii
    vdw_radii_X1 = vdw_radii[contact_indices[:, 0]]
    vdw_radii_X2 = vdw_radii[contact_indices[:, 1]]
    vdw_sum = (vdw_radii_X1 + vdw_radii_X2) / 10

    dists = calc_dist(xyz[contact_indices[:, 0]], xyz[contact_indices[:, 1]])
    vdw_contact = dists <= vdw_sum
    contact_indices = contact_indices[vdw_contact]

    return contact_indices


@njit(parallel=False, fastmath=True)
def dha_contacts(xyz, k, contact_id,
                 s1_donors, s1_hydros, s2_acc,
                 cut_HA, cut_DA, cut_DHA):
    """
    Args:
        xyz: Coordinates of the atoms in the frame.
        k (int): Frame identifier.
        contact_id: Identifier for the contact type.

        s1_donors: Indices of the donor atoms in s1.
        s1_hydros: Indices of the hydrogen atoms in s1.
        s2_acc: Indices of the acceptor atoms in s2.

        cut_HA: Cutoff distance between Hydrogen and Acceptor.
        cut_DA: Cuttoff distance between Donor and Acceptor.
        cut_DHA: Cutoff angle between Donor, Hydrogen and Acceptor.

    Returns:
        DHA_labels (list): List of DHA labels for each hydrogen bond.
    """

    # Create & query the trees
    s2_A_tree = nckd(xyz[s2_acc])
    ball_1 = s2_A_tree.query_radius(xyz[s1_hydros], cut_HA)

    # Find the DHA candidates
    n_triples = 0
    for x in ball_1:
        n_triples += x.size
    DHA_labels = np.empty((n_triples, 4), dtype=np.int32)
    DHA_coords = np.empty((n_triples, 3, 3), dtype=np.float32)

    counter = 0
    for i in range(len(ball_1)):
        D = s1_donors[i]
        H = s1_hydros[i]
        for j in ball_1[i]:
            A = s2_acc[j]
            DHA_labels[counter] = [D, A, k, contact_id]
            DHA_coords[counter][0] = xyz[D]
            DHA_coords[counter][1] = xyz[H]
            DHA_coords[counter][2] = xyz[A]
            counter += 1

    # Filter the DHA candidates
    DA_dist_correct = calc_dist(DHA_coords[:, 0],
                                DHA_coords[:, 2]) <= cut_DA
    DHA_angle_correct = calc_angle(DHA_coords[:, 0, :],
                                   DHA_coords[:, 1, :],
                                   DHA_coords[:, 2, :]) > cut_DHA

    correct_DHA = DHA_labels[DA_dist_correct & DHA_angle_correct]

    return correct_DHA.astype(np.int32)


# %% ==========================================================================
# Instantiate the IndexManager class
# =============================================================================
# todo: implement pi-cation
# todo: keep an eye on selections for contacts calculation (directed vs undirected)
# todo: watch angular definitions for the interactions
# todo: assign pi-type based on the angle between the planes
# todo: test interactions
# todo: reorganize code
# todo: parallelize
# todo: benchmark

import time
from intermap.indices import IndexManager

start_time = time.time()
topo = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
traj = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
sel1 = "nucleic or resname 8OG"
sel2 = "protein"
inters = 'all'

iman = IndexManager(topo, traj, sel1, sel2, inters)
load_time = time.time() - start_time
print(f"Until setting indices via SMARTS: {load_time:.2f} s")

# =============================================================================
# Load coordinates
# =============================================================================
import intermap.topo_trajs as tt

# Trajectory load
start = 0
last = 50
stride = 1
chunk_size = (last - start) // 1
chunks = tt.get_traj_chunks(topo, [traj],
                            start=start, last=last, stride=stride,
                            chunk_size=chunk_size)
chunk = next(chunks)
xyz_all = chunk.xyz
xyz = chunk.xyz[0]
k = 0

chunk_time = time.time() - start_time
print(f"Until loading the chunk: {chunk_time:.2f} s")


# =============================================================================
# unified approach
# =============================================================================


@njit(parallel=False)
def indices(full, subset):
    """
    Find the indices of the subset elements in the full array.

    Args:
        full (ndarray): Full array.
        subset (ndarray): Subset array.

    Returns:
        indices (ndarray): Indices of the subset elements in the full array.
    """
    indices = np.full(len(subset), -1, dtype=np.int32)
    for i in range(len(subset)):
        for j in range(len(full)):
            if full[j] == subset[i]:
                indices[i] = j
                break
    return indices


@njit(parallel=False)
def isin(a, b):
    n = len(a)
    result = np.full(n, False)
    set_b = set(b)
    for i in range(n):
        if a[i] in set_b:
            result[i] = True
    return result


@njit(parallel=False)
def duplex(xyz, k, inter_ids, s1_indices, s2_indices, dist_cut, anions,
           cations, hydroph, metal_don, metal_acc, vdw_radii, hb_hydros,
           hb_don, hb_acc):
    # Create & query the trees
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], dist_cut)

    # Find the close contacts
    n_contacts = sum([len(x) for x in ball_1])
    ijf = np.zeros((n_contacts, 3), dtype=np.int32)

    # Fill the ijf array
    counter = -1
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            counter += 1
            ijf[counter][0] = X1  # i
            ijf[counter][1] = X2  # j
            ijf[counter][2] = k  # f

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = ijf[:, 0] == ijf[:, 1]
    if idems.any():
        ijf = ijf[~idems]

    # Compute distances & allocate angles
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    n_types = len(inter_ids)
    interactions = np.zeros((n_contacts, n_types), dtype=np.bool_)

    # [CLOSE CONTACTS]
    cc_cut = cf.close_contacts
    cc_dist = dists <= cc_cut
    interactions[:, inter_ids['close_contacts']] = cc_dist

    # [ANIONIC]
    ionic_cut = cf.ionic
    ionic_dist = dists <= ionic_cut

    s1_ani = isin(row1, anions)
    s2_cat = isin(row2, cations)
    ani = s1_ani & s2_cat
    interactions[:, inter_ids['anionic']] = ani & ionic_dist

    # [CATIONIC]
    s1_cat = isin(row1, cations)
    s2_ani = isin(row2, anions)
    cat = s1_cat & s2_ani
    interactions[:, inter_ids['cationic']] = cat & ionic_dist

    # [HYDROPHOBIC]
    hp_cut = cf.hydrophobic
    hp_dist = dists <= hp_cut

    s1_hp = isin(row1, hydroph)
    s2_hp = isin(row2, hydroph)
    hp = s1_hp & s2_hp
    interactions[:, inter_ids['hydrophobic']] = hp & hp_dist

    # [METAL DONOR]
    met_cut = cf.metallic
    met_dist = dists <= met_cut

    s1_met_donors = isin(row1, metal_don)
    s2_met_acc = isin(row2, metal_acc)
    met_donors = s1_met_donors & s2_met_acc
    interactions[:, inter_ids['metal_donor']] = met_donors & met_dist

    # [METAL ACCEPTOR]
    s1_met_acc = isin(row1, metal_acc)
    s2_met_don = isin(row2, metal_don)
    met_acc = s1_met_acc & s2_met_don
    interactions[:, inter_ids['metal_acceptor']] = met_acc & met_dist

    # [VDW CONTACTS]
    row1_vdw = vdw_radii[row1]
    row2_vdw = vdw_radii[row2]
    vdw_sum = (row1_vdw + row2_vdw) / 10
    vdw_dist = dists <= vdw_sum
    interactions[:, inter_ids['vdw_contacts']] = vdw_dist

    # [HBOND ACCEPTOR]
    hb_ha_cut = cf.HA_cut
    hb_da_cut = cf.DA_cut
    hb_dha_cut = cf.DHA_cut
    hb_ha_dist = dists <= hb_ha_cut
    hb_da_dist = dists <= hb_da_cut
    hb_dists = hb_ha_dist & hb_da_dist

    s1_hb_acc = isin(row1, hb_acc)
    s2_hb_hydros = isin(row2, hb_hydros)
    before_angle = hb_dists & s1_hb_acc & s2_hb_hydros

    hydros = row2[before_angle]
    acc = row1[before_angle]
    idx = indices(hb_hydros, hydros)
    donors = hb_donors[idx]
    angles = calc_angle(xyz[donors], xyz[hydros], xyz[acc])

    count = 0
    for i in range(before_angle.size):
        if before_angle[i]:
            interactions[i, inter_ids['hb_acceptor']] = angles[
                                                            count] > hb_dha_cut
            count += 1

    # [HBOND DONORS]
    s1_hb_hydros = isin(row1, hb_hydros)
    s2_hb_acc = isin(row2, hb_acc)
    before_angle = hb_dists & s1_hb_hydros & s2_hb_acc

    hydros = row1[before_angle]
    acc = row2[before_angle]
    idx = indices(hb_acc, acc)
    donors = hb_don[idx]
    angles = calc_angle(xyz[donors], xyz[hydros], xyz[acc])

    count = 0
    for i in range(before_angle.size):
        if before_angle[i]:
            interactions[i, inter_ids['hb_donor']] = angles[count] > hb_dha_cut
            count += 1
    return ijf, dists, interactions


max_vdw = iman.get_max_vdw_dist()
dist_cut = max([
    cf.close_contacts,
    cf.ionic,
    cf.hydrophobic,
    cf.metallic,
    cf.hydrophobic,
    max_vdw])

s1_indices = iman.sel1_idx
s2_indices = iman.sel2_idx
anions = iman.anions
cations = iman.cations
hydroph = iman.hydroph
metal_don = iman.metal_don
metal_acc = iman.metal_acc
inter_ids = create_inter_ids()
vdw_radii = iman.radii
hb_hydros = iman.hb_H
hb_donors = iman.hb_D
hb_acc = iman.hb_A

ijf, dists, inters = duplex(xyz, k, inter_ids, s1_indices, s2_indices,
                            dist_cut, anions, cations, hydroph, metal_don,
                            metal_acc, vdw_radii, hb_hydros, hb_donors, hb_acc)

# =============================================================================
# %% Start computing interactions
# =============================================================================
start = time.time()
import intermap.cutoffs as cf

# ==== close contacts =========================================================
s1_indices = iman.sel1_idx
s2_indices = iman.sel2_idx
dist_cut = cf.close_contacts
cc_id = inter_identifiers['close_contacts']
if s1_indices.size and s2_indices.size:
    cc = close_contacts(xyz, k, cc_id, s1_indices, s2_indices, dist_cut)

# ==== anionic ================================================================
s1_ani = np.intersect1d(iman.sel1_idx, iman.anions)
s2_cat = np.intersect1d(iman.sel2_idx, iman.cations)
anionic_id = inter_identifiers['anionic']
if s1_ani.size and s2_cat.size:
    ani = close_contacts(xyz, k, anionic_id, s1_ani, s2_cat, cf.ionic)

# ==== cationic ===============================================================
s1_cat = np.intersect1d(iman.sel1_idx, iman.cations)
s2_ani = np.intersect1d(iman.sel2_idx, iman.anions)
cationic_id = inter_identifiers['cationic']
if s1_cat.size and s2_ani.size:
    cat = close_contacts(xyz, k, cationic_id, s1_cat, s2_ani, cf.ionic)

# ==== hydrophobic ============================================================
s1_hp = np.intersect1d(iman.sel1_idx, iman.hydroph)
s2_hp = np.intersect1d(iman.sel2_idx, iman.hydroph)
hydrop_id = inter_identifiers['hydrophobic']
if s1_hp.size and s2_hp.size:
    hp = close_contacts(xyz, k, hydrop_id, s1_hp, s2_hp, cf.hydrophobic)

# ==== metal donor ============================================================
s1_met_donors = np.intersect1d(iman.sel1_idx, iman.metal_don)
s2_met_acc = np.intersect1d(iman.sel2_idx, iman.metal_acc)
metd_id = inter_identifiers['metal_donor']
if s1_met_donors.size and s2_met_acc.size:
    metd = close_contacts(xyz, k, metd_id, s1_met_donors, s2_met_acc,
                          cf.metallic)

# ==== metal acceptor =========================================================
s1_acc = np.intersect1d(iman.sel1_idx, iman.metal_acc)
s2_met = np.intersect1d(iman.sel2_idx, iman.metal_don)
meta_id = inter_identifiers['metal_acceptor']
if s1_acc.size and s2_met.size:
    meta = close_contacts(xyz, k, meta_id, s1_acc, s2_met, cf.metallic)

# ==== vdw contacts ===========================================================
max_vdw = iman.get_max_vdw_dist()
vdw_radii = iman.radii
vdw_id = inter_identifiers['vdw_contacts']
if s1_indices.size and s2_indices.size:
    vdw = close_contacts(xyz, k, vdw_id, s1_indices, s2_indices, max_vdw,
                         vdw_radii)

# ==== hbonds acceptor ========================================================
idx = np.isin(iman.hb_H, iman.sel2_idx)
s2_hb_hydros = iman.hb_H[idx]
s2_hb_donors = iman.hb_D[idx]
s1_hb_acc = iman.hb_A[np.isin(iman.hb_A, iman.sel1_idx)]
hb_acc_id = inter_identifiers['hb_acceptor']
if s2_hb_donors.size and s1_hb_acc.size:
    hb_a = dha_contacts(xyz, k, hb_acc_id,
                        s2_hb_donors, s2_hb_hydros, s1_hb_acc,
                        cf.HA_cut, cf.DA_cut, cf.DHA_cut)

# ==== hbonds donor ===========================================================
idx = np.isin(iman.hb_H, iman.sel1_idx)
s1_hb_hydros = iman.hb_H[idx]
s1_hb_donors = iman.hb_D[idx]
s2_hb_acc = iman.hb_A[np.isin(iman.hb_A, iman.sel2_idx)]
hb_don_id = inter_identifiers['hb_donor']
if s1_hb_donors.size and s2_hb_acc.size:
    hb_d = dha_contacts(xyz, k, hb_don_id,
                        s1_hb_donors, s1_hb_hydros, s2_hb_acc,
                        cf.HA_cut, cf.DA_cut, cf.DHA_cut)
print(f"Computing time: {time.time() - start:.2f} s")

# if s1_donors.size and s1_hydros.size and s1_acc.size and \
#         s2_donors.size and s2_hydros.size and s2_acc.size:
#     hb = fnd.double_contacts(xyz, k, 0,
#                              s1_donors, s1_hydros, s1_acc,
#                              s2_donors, s2_hydros, s2_acc,
#                              cf.HA_cut, cf.DA_cut, cf.DHA_cut)
#
# if self.xb_D.size and self.xb_H.size and self.xb_A.size:
#     xb = fnd.double_contacts(xyz, k, 1,
#                              self.xb_D, self.xb_H, self.xb_A,
#                              self.xb_D, self.xb_H, self.xb_A,
#                              cf.HA_cut, cf.DA_cut, cf.DHA_cut)
#
# if sel1_rings.size and sel2_rings.size:
#     pipi = fnd.pi_stacking(xyz, k, 7, sel1_rings, sel2_rings, 0.6, 0.38, 0, 90)
#
# print(f"Computing time: {time.time() - start:.2f} s")


# =============================================================================
# Parallelize
# =============================================================================
import numpy as np
from numba import njit, prange


@njit(parallel=True)
def get_estimation(xyz_all, n_samples, factor=2):
    # Preallocate the arrays
    orig_len = xyz_all.shape[0]
    xyz_samples = xyz_all[::orig_len // n_samples]
    real_counts = np.zeros(n_samples, dtype=np.int32)

    # Use parallel loop to fill
    for i in prange(xyz_samples.shape[0]):
        cc = close_contacts(xyz_samples[i], 0, 0,
                            s1_indices, s2_indices,
                            dist_cut)
        real_counts[i] = cc.shape[0]

    # Estimate the number of contacts in the processed chunk
    estimation = np.int32(real_counts.max() * factor)
    return estimation


@njit(parallel=True)
def testing(xyz_all, n_samples):
    # Estimate the number of contacts in the processed chunk
    estimated = get_estimation(xyz_all, n_samples)
    # Preallocate the arrays
    num_frames = xyz_all.shape[0]
    all_cc = np.zeros((num_frames, estimated, 4), dtype=np.int32)
    real_counts = np.zeros(num_frames, dtype=np.int32)

    # Use parallel loop to fill
    for i in prange(xyz_all.shape[0]):
        cc = close_contacts(xyz_all[i], 0, 0, s1_indices, s2_indices, dist_cut)
        all_cc[i, :cc.shape[0]] = cc
        real_counts[i] = cc.shape[0]

    # Select the real counts
    all_cc_real = np.zeros((real_counts.sum(), 4), dtype=np.int32)
    counter = 0
    for i in range(num_frames):
        n = real_counts[i]
        all_cc_real[counter:counter + n] = all_cc[i, :n]
        counter += n
    return real_counts, all_cc, all_cc_real


real_counts, empty, good = testing(xyz_all, 5)
