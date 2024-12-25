# Created by rglez at 12/10/24

import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.commons as cmn
import intermap.cutoffs as cf


@njit(parallel=False)
def detect_vdw(inter_name, vdw_radii, row1, row2, dists):
    """
    Detect the Van der Waals interactions

    Args:
        inter_name (str): Name of the interaction
        vdw_radii (ndarray): Van der Waals radii of the atoms
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        dists (ndarray): Distances between the atoms in the first and second
                         selections

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        passing_dist (ndarray): Boolean array with the atoms that pass the
                                distance criterion
    """
    inter_idx = cmn.indices(to_compute_others, [inter_name])[0]

    if row1.size == 0:
        return inter_idx, np.zeros(dists.size, dtype=np.bool_)

    row1_vdw = vdw_radii[row1]
    row2_vdw = vdw_radii[row2]
    vdw_sum = (row1_vdw + row2_vdw) / 10
    passing_dist = dists <= vdw_sum
    return inter_idx, passing_dist


@njit(parallel=False)
def detect_1d(inter_name, dists, row1, type1, row2, type2):
    """
    Detect the 1D interactions (only one distance involved)

    Args:
        inter_name (str): Name of the interaction
        dists (ndarray): Distances between the atoms in the first and second
        row1 (ndarray): Indices of the atoms in the first selection
        type1 (ndarray): Indices of the atoms of the first type
        row2 (ndarray): Indices of the atoms in the second selection
        type2 (ndarray): Indices of the atoms of the second type

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        passing_dist (ndarray): Boolean array with the atoms that pass the
                                distance criterion
    """
    inter_idx = cmn.indices(to_compute_others, [inter_name])[0]

    if type1.size == 0 or type2.size == 0:
        return inter_idx, np.zeros(dists.size, dtype=np.bool_)

    dist_cutoff = cutoffs_others[0, inter_idx]
    passing_dists = dists <= dist_cutoff

    if inter_name == 'CloseContacts':
        return inter_idx, passing_dists

    else:
        s1_is_type = cmn.isin(row1, type1)
        s2_is_type = cmn.isin(row2, type2)
        are_type = s1_is_type & s2_is_type
        return inter_idx, passing_dists & are_type


@njit(parallel=False)
def detect_2d1a(inter_name, dists, row1, row2, hb_acc, hb_hydros, hb_donors,
                ha_cut, da_cut, min_ang, max_ang):
    """
    Detect the 2D1A interactions (two distances and one angle involved)

    Args:
        inter_name (str): Name of the interaction
        dists (ndarray): Distances between the atoms in the first and second
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        hb_acc (ndarray): Indices of the hydrogen bond acceptors
        hb_hydros (ndarray): Indices of the hydrogen bond hydrogens
        hb_donors (ndarray): Indices of the hydrogen bond donors
        ha_cut (float): Cutoff distance for the hydrogen bond acceptor
        da_cut (float): Cutoff distance for the hydrogen bond donor-acceptor
        min_ang (float): Minimum angle for the hydrogen bond
        max_ang (float): Maximum angle for the hydrogen bond

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        hbonds (ndarray): Boolean array with the atoms that pass the distance
    """
    idx_name = cmn.indices(to_compute_others, [inter_name])[0]

    if 'Acceptor' in inter_name:
        are_acceptors = cmn.isin(row1, hb_acc)
        are_hydros = cmn.isin(row2, hb_hydros)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        acceptors = row1[passing_ha]
        hydros = row2[passing_ha]

    elif 'Donor' in inter_name:
        are_hydros = cmn.isin(row1, hb_hydros)
        are_acceptors = cmn.isin(row2, hb_acc)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        hydros = row1[passing_ha]
        acceptors = row2[passing_ha]
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")

    if passing_ha.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    ext_idx = cmn.indices(hb_hydros, hydros)
    donors = hb_donors[ext_idx]

    ha_vectors = xyz[acceptors] - xyz[hydros]
    hd_vectors = xyz[donors] - xyz[hydros]
    da_vectors = xyz[donors] - xyz[acceptors]

    da_dists = np.sqrt((da_vectors ** 2).sum(axis=1))
    passing_da = da_dists <= da_cut

    if not passing_da.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    angles = cmn.calc_angles_2v(hd_vectors[passing_da], ha_vectors[passing_da])
    passing_angles = (angles >= min_ang) & (angles <= max_ang)

    if passing_angles.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    hbonds = np.zeros(dists.size, dtype=np.bool_)
    hbonds[passing_ha] = passing_angles
    return idx_name, hbonds


@njit(parallel=False)
def stackings(inter_name, dists, mindists, s1_normals, s2_normals):
    """
    Helper function to compute the pi-stacking interactions

    """

    # Parse the cutoffs
    idx = cmn.indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_dist = cutoffs_aro[1, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Apply restraints
    passing_dist1 = dists <= dist_cut
    passing_dist2 = mindists <= min_dist
    angles = cmn.calc_angles_2v(s1_normals, s2_normals)
    passing_angles = ((angles >= min_ang) & (angles <= max_ang))
    stacking = passing_dist1 & passing_dist2 & passing_angles
    return idx, stacking


@njit(parallel=False)
def pications(inter_name, ijf, dists, xyz2, rings_normals, rings_idx, cat_idx,
              to_compute_aro):
    """
    Helper function to compute the pi-cation // cation-pi interactions

    """
    # Parse the cutoffs
    idx = cmn.indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Select the pairs
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    if inter_name == 'PiCation':
        s1_is_type = cmn.isin(row1, rings_idx)
        s2_is_type = cmn.isin(row2, cat_idx)
    elif inter_name == 'CationPi':
        s1_is_type = cmn.isin(row1, cat_idx)
        s2_is_type = cmn.isin(row2, rings_idx)
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")
    pairs = s1_is_type & s2_is_type

    if pairs.any():
        # Calculate angles between normals and vectors
        row1_pairs = row1[pairs]
        row2_pairs = row2[pairs]
        vector_ctr_cat = xyz2[row1_pairs] - xyz2[row2_pairs]
        normals = rings_normals[cmn.indices(rings_idx, row1_pairs)]
        angles = cmn.calc_angles_2v(normals, vector_ctr_cat)

        # Apply restraints
        passing_dist = dists[pairs] <= dist_cut
        passing_angles = ((angles >= min_ang) & (angles <= max_ang))
        return idx, pairs, passing_dist & passing_angles
    else:
        return idx, pairs, pairs


@njit(parallel=False)
def not_aro(xyz, k, s1_indices_raw, s2_indices_raw, anions,
            cations, hydroph, metal_don, metal_acc, vdw_radii, hb_hydros,
            hb_donors, hb_acc, cutoffs_others, to_compute_others):
    """
    Compute the interactions not related to aromatic rings

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices_raw (ndarray): Indices of the atoms in the first selection
        s2_indices_raw (ndarray): Indices of the atoms in the second selection
        anions (ndarray): Indices of the anions
        cations (ndarray): Indices of the cations
        hydroph (ndarray): Indices of the hydrophobic atoms
        metal_don (ndarray): Indices of the metal donors
        metal_acc (ndarray): Indices of the metal acceptors
        vdw_radii (ndarray): Van der Waals radii of the atoms
        hb_hydros (ndarray): Indices of the hydrogen bond hydrogens
        hb_donors (ndarray): Indices of the hydrogen bond donors
        hb_acc (ndarray): Indices of the hydrogen bond acceptors
        cutoffs_others (ndarray): Cutoff distances for the interactions not
        to_compute_others (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    # Create & query the trees
    dist_cut = cutoffs_others[:2].max()
    s2_tree = nckd(xyz[s2_indices_raw])
    ball_1 = s2_tree.query_radius_parallel(xyz[s1_indices_raw], dist_cut, p=np.inf)

    # Get containers
    ijf, dists, interactions = cmn.get_containers(
        xyz, k, np.arange(xyz.shape[0]), ball_1, s1_indices_raw,
        s2_indices_raw, to_compute_others)

    # Start detecting interactions
    to_detect = list(to_compute_others)
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]

    # Contacts
    inter_name = 'VdWContact'
    if inter_name in to_detect:
        inter_idx, vdw = detect_vdw(inter_name, vdw_radii, row1, row2, dists)
        interactions[:, inter_idx] = vdw

    inter_name = 'CloseContacts'
    if inter_name in to_detect:
        inter_idx, close_contacts = detect_1d(inter_name, dists, row1, row1,
                                              row2, row2)
        interactions[:, inter_idx] = close_contacts

    inter_name = 'Hydrophobic'
    if inter_name in to_detect:
        inter_idx, hydrophobic = detect_1d(inter_name, dists, row1, hydroph,
                                           row2, hydroph)
        interactions[:, inter_idx] = hydrophobic

    # Salt bridges
    exist_all = anions.size > 0 and cations.size > 0
    find_cationic = ('Cationic' in to_detect) and exist_all
    find_anionic = ('Anionic' in to_detect) and exist_all

    if find_cationic:
        inter_name = 'Cationic'
        inter_idx, cationic = detect_1d(inter_name, dists, row1, cations, row2,
                                        anions)
        interactions[:, inter_idx] = cationic

    if find_anionic in to_detect:
        inter_name = 'Anionic'
        inter_idx, anionic = detect_1d(inter_name, dists, row1, anions, row2,
                                       cations)
        interactions[:, inter_idx] = anionic

    # Metallic
    exist_all = metal_don.size > 0 and metal_acc.size > 0
    find_metal_don = ('MetalDonor' in to_detect) and exist_all
    find_metal_acc = ('MetalAcceptor' in to_detect) and exist_all

    if find_metal_don in to_detect:
        inter_name = 'MetalDonor'
        inter_idx, metal_donor = detect_1d(inter_name, dists, row1, metal_don,
                                           row2, metal_acc)
        interactions[:, inter_idx] = metal_donor

    if find_metal_acc in to_detect:
        inter_name = 'MetalAcceptor'
        inter_idx, metal_acc = detect_1d(inter_name, dists, row1, metal_acc,
                                         row2, metal_don)
        interactions[:, inter_idx] = metal_acc

    # HBonds
    exists_all = hb_acc.size > 0 and hb_donors.size > 0
    find_hb_acc = ('HBAcceptor' in to_detect) and exists_all
    find_hb_donor = ('HBDonor' in to_detect) and exists_all
    if find_hb_acc or find_hb_donor:
        idx_both = cmn.indices(to_compute, ['HBAcceptor'])[0]
        da_cut, ha_cut, min_ang, max_ang = cutoffs_others[:4, idx_both]

        if find_hb_acc:
            inter_name = 'HBAcceptor'
            inter_idx, hb_a = detect_2d1a(inter_name, dists, row1, row2,
                                          hb_acc, hb_hydros, hb_donors,
                                          ha_cut, da_cut, min_ang, max_ang)
            interactions[:, inter_idx] = hb_a

        if find_hb_donor:
            inter_name = 'HBDonor'
            inter_idx, hb_d = detect_2d1a(inter_name, dists, row1, row2,
                                          hb_acc, hb_hydros, hb_donors, ha_cut,
                                          da_cut, min_ang, max_ang)
            interactions[:, inter_idx] = hb_d

    # XBonds
    exists_all = xb_acceptors.size > 0 and xb_donors.size > 0
    find_xb_acc = ('XBAcceptor' in to_detect) and exists_all
    find_xb_donor = ('XBDonor' in to_detect) and exists_all

    if find_xb_acc or find_xb_donor:
        idx_both = cmn.indices(to_compute, ['XBAcceptor'])[0]
        da_cut, ha_cut, min_ang, max_ang = cutoffs_others[:4, idx_both]

        if find_xb_acc:
            inter_name = 'XBAcceptor'
            inter_idx, xb_a = detect_2d1a(inter_name, dists, row1, row2,
                                          xb_acceptors, xb_halogens, xb_donors,
                                          ha_cut, da_cut, min_ang, max_ang)
            interactions[:, inter_idx] = xb_a

        if find_xb_donor:
            inter_name = 'XBDonor'
            inter_idx, xb_d = detect_2d1a(inter_name, dists, row1, row2,
                                          xb_acceptors, xb_halogens, xb_donors,
                                          ha_cut, da_cut, min_ang, max_ang)
            interactions[:, inter_idx] = xb_d

    return ijf, interactions


@njit(parallel=False)
def aro(xyz, k, s1_indices_raw, s2_indices_raw, cations, rings, cutoffs_aro,
        to_compute_aro):
    """
    Compute the aromatic interactions

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices_raw (ndarray): Indices of the atoms in the first selection
        s2_indices_raw (ndarray): Indices of the atoms in the second selection
        cations (ndarray): Indices of the cations
        rings (ndarray): Indices of the aromatic rings
        cutoffs_aro (ndarray): Cutoff distances for the aromatic interactions
        to_compute_aro (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """
    # =========================================================================
    # STEP I: Find all pair of atoms/centroids within the max cutoff distance
    # =========================================================================

    # Get cations
    s1_cat = s1_indices_raw[cmn.isin(s1_indices_raw, cations)]
    s2_cat = s2_indices_raw[cmn.isin(s2_indices_raw, cations)]

    # Get the aromatic rings
    s1_rings = rings[cmn.isin(rings[:, 0], s1_indices_raw)]
    s2_rings = rings[cmn.isin(rings[:, 0], s2_indices_raw)]

    # Compute the centroids
    s1_centr = cmn.calc_centroids(s1_rings, xyz)
    s2_centr = cmn.calc_centroids(s2_rings, xyz)

    # Compute the normal vectors
    s1_at1, s1_at3, s1_at5 = s1_rings[:, 0], s1_rings[:, 2], s1_rings[:, 4]
    s2_at1, s2_at3, s2_at5 = s2_rings[:, 0], s2_rings[:, 2], s2_rings[:, 4]
    s1_norm = cmn.calc_normal_vector(xyz[s1_at1], xyz[s1_at3], xyz[s1_at5])
    s2_norm = cmn.calc_normal_vector(xyz[s2_at1], xyz[s2_at3], xyz[s2_at5])

    # Create a new xyz array with the cations & centroids only
    xyz2 = concat((xyz[s1_cat], xyz[s2_cat], s1_centr, s2_centr), axis=0)
    ext_idx = concat((s1_cat, s2_cat, s1_rings[:, 0], s2_rings[:, 0]))

    # Internal indexing for xyz2 coordinates
    n0 = s1_cat.size + s2_cat.size
    n1 = n0 + s1_centr.shape[0]
    n2 = n1 + s2_centr.shape[0]
    s1_cat_idx = np.arange(0, s1_cat.size)
    s2_cat_idx = np.arange(s1_cat.size, n0)
    s1_rings_idx = np.arange(n0, n1)
    s2_rings_idx = np.arange(n1, n2)

    # Create & query the trees
    s1_indices = concat((s1_cat_idx, s1_rings_idx))
    s2_indices = concat((s2_cat_idx, s2_rings_idx))
    dist_cut_aro = cutoffs_aro[:2].max()

    if s2_indices.size == 0 or s1_indices.size == 0:
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, 1), dtype=np.bool_)
        return ijf, inters

    s2_tree = nckd(xyz2[s2_indices])
    ball_1 = s2_tree.query_radius_parallel(xyz2[s1_indices], dist_cut_aro, p=2)
    ijf, dists, inters = cmn.get_containers(
        xyz2, k, ext_idx, ball_1, s1_indices, s2_indices, to_compute_aro)

    # =========================================================================
    # STEP II: Compute the aromatic interactions
    # =========================================================================
    set_aro = list(to_compute_aro)

    if 'PiCation' in set_aro:
        idx, pairs, pi_cat = pications('PiCation', ijf, dists, xyz2, s1_norm,
                                       s1_rings_idx, s2_cat_idx,
                                       to_compute_aro)
        if pairs.any():
            inters[pairs, idx] = pi_cat

    if 'CationPi' in set_aro:
        idx, pairs, cat_pi = pications('CationPi', ijf, dists, xyz2, s2_norm,
                                       s2_rings_idx, s1_cat_idx,
                                       to_compute_aro)
        if pairs.any():
            inters[pairs, idx] = cat_pi

    # Stackings
    find_PiStacking = 'PiStacking' in set_aro
    find_EdgeToFace = 'EdgeToFace' in set_aro
    find_FaceToFace = 'FaceToFace' in set_aro
    if find_PiStacking or find_EdgeToFace or find_FaceToFace:

        # Get the ring pairs
        row1, row2 = ijf[:, 0], ijf[:, 1]
        s1_is_ctr = cmn.isin(row1, s1_rings_idx)
        s2_is_ctr = cmn.isin(row2, s2_rings_idx)
        pairs = s1_is_ctr & s2_is_ctr
        if pairs.any():
            ring_pairs = ijf[pairs]
            ring_dists = dists[pairs]
            s1_normals = s1_norm[cmn.indices(s1_rings_idx, ring_pairs[:, 0])]
            s2_normals = s2_norm[cmn.indices(s2_rings_idx, ring_pairs[:, 1])]

            # Compute the minimum distance between the rings
            num_pairs = ring_pairs.shape[0]
            mindists = np.zeros(num_pairs, dtype=np.float32)
            for i in range(num_pairs):
                s1_ring = \
                    s1_rings[cmn.indices(s1_rings_idx, [ring_pairs[i, 0]])][0]
                s1_ring_idx = s1_ring[:s1_ring[-1]]

                s2_ring = \
                    s2_rings[cmn.indices(s2_rings_idx, [ring_pairs[i, 1]])][0]
                s2_ring_idx = s2_ring[:s2_ring[-1]]
                mindists[i] = cmn.calc_min_dist(xyz[s1_ring_idx],
                                                xyz[s2_ring_idx])

            if find_PiStacking:
                idx, pi_stacking = stackings('PiStacking', ring_dists,
                                             mindists,
                                             s1_normals, s2_normals)
                inters[pairs, idx] = pi_stacking

            if find_EdgeToFace:
                idx, etf_stacking = stackings('EdgeToFace', ring_dists,
                                              mindists,
                                              s1_normals, s2_normals)
                inters[pairs, idx] = etf_stacking

            if find_FaceToFace:
                idx, ftf_stacking = stackings('FaceToFace', ring_dists,
                                              mindists,
                                              s1_normals, s2_normals)
                inters[pairs, idx] = ftf_stacking
    return ijf, inters


# %% ==========================================================================
# Debugging Area
# =============================================================================

import time
from intermap.indices import IndexManager

start_time = time.time()
topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
sel1 = "nucleic"
sel2 = "resname WAT"
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

# %% ==========================================================================
# Runner (to put there later)
# =============================================================================
import re

# Parse the cutoffs & select the interactions to compute
parsed_cutoffs = cf.parse_cutoffs(args=())
all_inters, all_cutoffs = cf.get_inters_cutoffs(parsed_cutoffs)

to_compute = 'all'
if to_compute == 'all':
    to_compute = all_inters
else:
    to_compute = to_compute

bit_aro = [i for i, x in enumerate(to_compute) if re.search(r'Pi|Face', x)]
bit_others = [i for i in range(all_cutoffs.shape[1]) if i not in bit_aro]
to_compute_aro = to_compute[bit_aro]
to_compute_others = to_compute[bit_others]
cutoffs_aro = all_cutoffs[:, bit_aro]
cutoffs_others = all_cutoffs[:, bit_others]

# Get the indices from the IndexManager object
s1_indices_raw = iman.sel1_idx
s2_indices_raw = iman.sel2_idx
vdw_radii = iman.radii
hydrophobes = iman.hydroph
anions = iman.anions
cations = iman.cations
metal_donors = iman.metal_don
metal_acceptors = iman.metal_acc
hb_hydrogens = iman.hb_H
hb_donors = iman.hb_D
hb_acceptors = iman.hb_A
xb_halogens = iman.xb_H
xb_donors = iman.xb_D
xb_acceptors = iman.xb_A
rings = iman.rings

# %% ==========================================================================
ijf_aro, inters_aro = aro(xyz, k, s1_indices_raw, s2_indices_raw, cations,rings, cutoffs_aro, to_compute_aro)
ijf_others, inters_others = not_aro(xyz, k, s1_indices_raw, s2_indices_raw,anions, cations, hydrophobes, metal_donors,metal_acceptors, vdw_radii, hb_hydrogens,hb_donors, hb_acceptors, cutoffs_others,to_compute_others)

%timeit ijf_aro, inters_aro = aro(xyz, k, s1_indices_raw, s2_indices_raw, cations,rings, cutoffs_aro, to_compute_aro)
%timeit ijf_others, inters_others = not_aro(xyz, k, s1_indices_raw, s2_indices_raw,anions, cations, hydrophobes, metal_donors,metal_acceptors, vdw_radii, hb_hydrogens,hb_donors, hb_acceptors, cutoffs_others,to_compute_others)

print('Aro', inters_aro.sum(axis=0))
print('Not-Aro', inters_others.sum(axis=0))

# =============================================================================
# Parallelize
# =============================================================================
# import numpy as np
# from numba import njit, prange
#
#
# @njit(parallel=True)
# def get_estimation(xyz_all, n_samples, factor=2):
#     # Preallocate the arrays
#     orig_len = xyz_all.shape[0]
#     xyz_samples = xyz_all[::orig_len // n_samples]
#     real_counts = np.zeros(n_samples, dtype=np.int32)
#
#     # Use parallel loop to fill
#     for i in prange(xyz_samples.shape[0]):
#         cc = close_contacts(xyz_samples[i], 0, 0,
#                             s1_indices, s2_indices,
#                             dist_cut)
#         real_counts[i] = cc.shape[0]
#
#     # Estimate the number of contacts in the processed chunk
#     estimation = np.int32(real_counts.max() * factor)
#     return estimation
#
#
# @njit(parallel=True)
# def testing(xyz_all, n_samples):
#     # Estimate the number of contacts in the processed chunk
#     estimated = get_estimation(xyz_all, n_samples)
#     # Preallocate the arrays
#     num_frames = xyz_all.shape[0]
#     all_cc = np.zeros((num_frames, estimated, 4), dtype=np.int32)
#     real_counts = np.zeros(num_frames, dtype=np.int32)
#
#     # Use parallel loop to fill
#     for i in prange(xyz_all.shape[0]):
#         cc = close_contacts(xyz_all[i], 0, 0, s1_indices, s2_indices, dist_cut)
#         all_cc[i, :cc.shape[0]] = cc
#         real_counts[i] = cc.shape[0]
#
#     # Select the real counts
#     all_cc_real = np.zeros((real_counts.sum(), 4), dtype=np.int32)
#     counter = 0
#     for i in range(num_frames):
#         n = real_counts[i]
#         all_cc_real[counter:counter + n] = all_cc[i, :n]
#         counter += n
#     return real_counts, all_cc, all_cc_real
#
#
# real_counts, empty, good = testing(xyz_all, 5)
