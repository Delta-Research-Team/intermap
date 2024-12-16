# Created by rglez at 12/10/24

import numba.types
import numpy as np
from numba import njit
from numba.typed import Dict as numba_dict
from numba_kdtree import KDTree as nckd

import intermap.cutoffs as cf


def parse_cutoff(cutoff_name, args):
    """
    Parse the cutoff name from the args or the cf module if not found in args.

    Args:
        args (Namespace): Arguments from the parser.
        cutoff_name (str): Name of the cutoff.

    Returns:
        float: Value of the cutoff.
    """
    if cutoff_name not in cf.__dict__:
        raise ValueError(f"{cutoff_name} is not a valid cutoff name.\n"
                         f"The supported list is:\n"
                         f"{[x for x in dir(cf) if not x.startswith("__")]}")
    elif cutoff_name in args:
        return getattr(args, cutoff_name)
    else:
        # get the value from the cf module
        return getattr(cf, cutoff_name)


@njit
def create_inter_ids():
    # Create a Numba Dict with specific types
    inter_ids = numba_dict.empty(numba.types.unicode_type, numba.types.int64)

    # undirected 2p interactions
    inter_ids['CloseContact'] = 0
    inter_ids['VdWContact'] = 1
    inter_ids['Hydrophobic'] = 2

    # directed 2p interactions
    inter_ids['Anionic'] = 3
    inter_ids['Cationic'] = 4
    inter_ids['MetalDonor'] = 5
    inter_ids['MetalAcceptor'] = 6

    # directed 3p interactions
    inter_ids['HBAcceptor'] = 7
    inter_ids['HBDonor'] = 8
    inter_ids['XBAcceptor'] = 9
    inter_ids['XBDonor'] = 10

    # undirected np interactions
    inter_ids['PiStacking'] = 11
    inter_ids['FaceToFace'] = 12
    inter_ids['EdgeToFace'] = 13

    # directed np interactions
    inter_ids['PiCation'] = 14
    inter_ids['CationPi'] = 15

    return inter_ids


interactions = {
    "Anionic": [parse_cutoff('dist_cut_Ionic')],
    "CationPi": [parse_cutoff('dist_cut_PiCation')],
    "Cationic": [parse_cutoff('dist_cut_Ionic')],
    "CloseContact": [parse_cutoff('dist_cut_CloseContacts')],
    "EdgeToFace": [parse_cutoff('dist_cut_EdgeToFace'),
                   parse_cutoff('min_dist_EdgeToFace')],
    "FaceToFace": [parse_cutoff('dist_cut_FaceToFace'),
                   parse_cutoff('min_dist_FaceToFace')],
    "HBAcceptor": [parse_cutoff('dist_cut_DA'), parse_cutoff('dist_cut_HA')],
    "HBDonor": [parse_cutoff('dist_cut_DA'), parse_cutoff('dist_cut_HA')],
    "Hydrophobic": [parse_cutoff('dist_cut_Hydroph')],
    "MetalAcceptor": [parse_cutoff('dist_cut_Metalic')],
    "MetalDonor": [parse_cutoff('dist_cut_Metalic')],
    "PiCation": [parse_cutoff('dist_cut_PiCation')],
    "PiStacking": [parse_cutoff('dist_cut_PiStacking'),
                   parse_cutoff('min_dist_PiStacking')],
    "VdWContact": [],
    "XBAcceptor": [parse_cutoff('dist_cut_XA'), parse_cutoff('dist_cut_XD')],
    "XBDonor": [parse_cutoff('dist_cut_XA'), parse_cutoff('dist_cut_XD')],
}


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
sel2 = sel1
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
def compute_centroids(rings, xyz):
    centroids = np.zeros((len(rings), 3), dtype=float)
    for i, ring in enumerate(rings):
        atoms = ring[:ring[-1]]
        centroids[i] = xyz[atoms].sum(axis=0) / len(atoms)
    return centroids


@njit(parallel=False)
def calc_normal_vector(p1, p2, p3):
    """
    Calculate the normal vector of a plane defined by three points

    Args:
        p1 (np.ndarray): First point
        p2 (np.ndarray): Second point
        p3 (np.ndarray): Third point

    Returns:
        np.ndarray: Normal vector of the plane defined by the three points
    """
    # Calculate vectors from points
    v1 = p2 - p1
    v2 = p3 - p1

    # Calculate the normal vector
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)

    return (normal / norm).astype(np.float32)


@njit(parallel=False)
def calc_min_dist(coords1, coords2):
    """
    Get the minimumm distance between two sets of coordinates

    Args:
        coords1: coordinates of the first residue
        coords2: coordinates of the second residue

    Returns:
        The minimum distance between two sets of coordinates
    """
    # Constants
    n1 = coords1.shape[0]
    n2 = coords2.shape[0]

    # Find minimum distance using square values to save time
    min_dist_squared = np.inf
    for i in range(n1):
        for j in range(n2):
            dist_squared = (
                    (coords1[i][0] - coords2[j][0]) ** 2
                    + (coords1[i][1] - coords2[j][1]) ** 2
                    + (coords1[i][2] - coords2[j][2]) ** 2
            )
            if dist_squared < min_dist_squared:
                min_dist_squared = dist_squared
    return np.sqrt(min_dist_squared)


@njit(parallel=False)
def duplex(xyz, k, inter_ids, s1_indices_raw, s2_indices_raw, dist_cut, anions,
           cations, hydroph, metal_don, metal_acc, vdw_radii, hb_hydros,
           hb_don, hb_acc, rings, to_compute):
    # =========================================================================
    # STEP I: Find all pair of atoms (and centroids) within the cutoff distance
    # =========================================================================

    # Add aromatic centroids to xyz
    s1_rings = rings[isin(rings[:, 0], s1_indices_raw)].astype(np.int32)
    s2_rings = rings[isin(rings[:, 0], s2_indices_raw)].astype(np.int32)
    sel1_centroids = compute_centroids(s1_rings, xyz)
    sel2_centroids = compute_centroids(s2_rings, xyz)
    xyz2 = np.concatenate((xyz, sel1_centroids, sel2_centroids), axis=0)

    # Update sel1 & sel2 indices with the aromatic centroids
    n0 = xyz.shape[0]
    n1 = n0 + sel1_centroids.shape[0]
    n2 = n1 + sel2_centroids.shape[0]
    s1_rings_indices = np.arange(n0, n1)
    s2_rings_indices = np.arange(n1, n2)
    s1_indices = np.concatenate((s1_indices_raw, s1_rings_indices))
    s2_indices = np.concatenate((s2_indices_raw, s2_rings_indices))

    # Create & query the trees
    s2_tree = nckd(xyz2[s2_indices])
    ball_1 = s2_tree.query_radius(xyz2[s1_indices], dist_cut)

    # Find the close contacts
    n_contacts = sum([len(x) for x in ball_1])
    ijf = np.zeros((n_contacts, 3), dtype=np.int32)
    counter = 0
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            ijf[counter][0] = X1
            ijf[counter][1] = X2
            ijf[counter][2] = k
            counter += 1

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = ijf[:, 0] == ijf[:, 1]
    if idems.any():
        ijf = ijf[~idems]

    # Compute distances
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = calc_dist(xyz2[row1], xyz2[row2])

    # Create the container for interaction types
    inter_ids_selected = [x for x in inter_ids.values() if x in to_compute]
    n_types = len(inter_ids_selected)
    interactions = np.zeros((ijf.shape[0], n_types), dtype=np.bool_)

    # =========================================================================
    # STEP II: Compute the interactions
    # =========================================================================

    # ---- [AROMATICS] --------------------------------------------------------
    sel1_normals = calc_normal_vector(xyz[s1_rings[:, 0]],
                                      xyz[s1_rings[:, 2]],
                                      xyz[s1_rings[:, 4]])
    sel2_normals = calc_normal_vector(xyz[s2_rings[:, 0]],
                                      xyz[s2_rings[:, 2]],
                                      xyz[s2_rings[:, 4]])

    # Cutoffs definitions
    stacking_cut = cf.dist_cut_PiStacking
    stacking_min = cf.min_ang_PiStacking
    stacking_max = cf.max_ang_PiStacking
    stacking_min_dist = cf.min_dist_PiStacking

    etf_cut = cf.dist_cut_EdgeToFace
    etf_min = cf.min_ang_EdgeToFace
    etf_max = cf.max_ang_EdgeToFace
    etf_min_dist = cf.min_dist_EdgeToFace

    ftf_cut = cf.dist_cut_FaceToFace
    ftf_min = cf.min_ang_FaceToFace
    ftf_max = cf.max_ang_FaceToFace
    ftf_min_dist = cf.min_dist_FaceToFace

    # Get the ring pairs
    s1_centroids = isin(row1, s1_rings_indices)
    s2_centroids = isin(row2, s2_rings_indices)
    close_centroids = s1_centroids & s2_centroids
    ring_pairs = ijf[close_centroids]

    # Compute the normal vectors
    s1_close_rings = s1_rings[indices(s1_rings_indices, ring_pairs[:, 0])]
    s2_close_rings = s2_rings[indices(s2_rings_indices, ring_pairs[:, 1])]
    s1_normals = calc_normal_vector(xyz2[s1_close_rings[:, 0]],
                                    xyz2[s1_close_rings[:, 2]],
                                    xyz2[s1_close_rings[:, 4]])
    s2_normals = calc_normal_vector(xyz2[s2_close_rings[:, 0]],
                                    xyz2[s2_close_rings[:, 2]],
                                    xyz2[s2_close_rings[:, 4]])

    # Compute the angle between the normals
    dot_products = np.einsum('ij,ij->i', s1_normals, s2_normals)
    dot_64 = np.float64(dot_products)
    dot_64[dot_64 < -1.0] = -1.0
    dot_64[dot_64 > 1.0] = 1.0
    radians = np.arccos(dot_64)
    angles = np.degrees(radians)

    # Compute the minimum distance between the rings
    num_pairs = ring_pairs.shape[0]
    mindists = np.zeros(num_pairs, dtype=np.float32)
    for i in range(num_pairs):
        r1 = s1_close_rings[i][:s1_close_rings[i][-1]]
        r2 = s2_close_rings[i][:s2_close_rings[i][-1]]
        mindists[i] = calc_min_dist(xyz2[r1], xyz2[r2])

    # pi-stacking interactions
    by_stacking_dist = dists[close_centroids] <= stacking_cut
    by_stacking_angles = (angles >= stacking_min) & (angles <= stacking_max)
    by_stacking_mindist = mindists <= stacking_min_dist
    pi_stacking = by_stacking_dist & by_stacking_angles & by_stacking_mindist
    interactions[close_centroids, inter_ids['pi_stacking']] = pi_stacking

    # edge-to-face interactions
    by_etf_dist = dists[close_centroids] <= etf_cut
    by_etf_angles = (angles >= etf_min) & (angles <= etf_max)
    by_etf_mindist = mindists <= etf_min_dist
    etf = by_etf_dist & by_etf_angles & by_etf_mindist
    interactions[close_centroids, inter_ids['edge_to_face']] = etf

    # face-to-face interactions
    by_ftf_dist = dists[close_centroids] <= ftf_cut
    by_ftf_angles = (angles >= ftf_min) & (angles <= ftf_max)
    by_ftf_mindist = mindists <= ftf_min_dist
    ftf = by_ftf_dist & by_ftf_angles & by_ftf_mindist
    interactions[close_centroids, inter_ids['face_to_face']] = ftf

    # ---- [Pi-Cation] --------------------------------------------------------
    # pi_cat_cut = cf.dist_cut_PiCation
    s2_cat = isin(row2, cations)
    # pi_cat = s1_centroids & s2_cat
    # pi_cat_dist = dists <= pi_cat_cut
    # compute_angles(s1)
    # interactions[:, inter_ids['pi_cation']] = pi_cat & pi_cat_dist

    # ---- [Cation-Pi] --------------------------------------------------------
    s1_cat = isin(row1, cations)
    # cat_pi = s1_cat & s2_centroids
    # cat_pi_dist = dists <= pi_cat_cut
    # interactions[:, inter_ids['cation_pi']] = cat_pi & cat_pi_dist

    # ---- [CLOSE CONTACTS] ---------------------------------------------------
    cc_cut = cf.dist_cut_CloseContacts
    cc_dist = dists <= cc_cut
    interactions[:, inter_ids['close_contacts']] = cc_dist

    # ---- [ANIONIC] ----------------------------------------------------------
    ionic_cut = cf.dist_cut_Ionic
    ionic_dist = dists <= ionic_cut

    s1_ani = isin(row1, anions)
    ani = s1_ani & s2_cat
    interactions[:, inter_ids['anionic']] = ani & ionic_dist

    # ---- [CATIONIC]- --------------------------------------------------------
    s2_ani = isin(row2, anions)
    cat = s1_cat & s2_ani
    interactions[:, inter_ids['cationic']] = cat & ionic_dist

    # ---- [HYDROPHOBIC] ------------------------------------------------------
    hp_cut = cf.dist_cut_Hydroph
    hp_dist = dists <= hp_cut

    s1_hp = isin(row1, hydroph)
    s2_hp = isin(row2, hydroph)
    hp = s1_hp & s2_hp
    interactions[:, inter_ids['hydrophobic']] = hp & hp_dist

    # ---- [METAL DONOR] ------------------------------------------------------
    met_cut = cf.dist_cut_Metalic
    met_dist = dists <= met_cut

    s1_met_donors = isin(row1, metal_don)
    s2_met_acc = isin(row2, metal_acc)
    met_donors = s1_met_donors & s2_met_acc
    interactions[:, inter_ids['metal_donor']] = met_donors & met_dist

    # ---- [METAL ACCEPTOR] ---------------------------------------------------
    s1_met_acc = isin(row1, metal_acc)
    s2_met_don = isin(row2, metal_don)
    met_acc = s1_met_acc & s2_met_don
    interactions[:, inter_ids['metal_acceptor']] = met_acc & met_dist

    # ---- [VDW CONTACTS] -----------------------------------------------------
    row1_vdw = vdw_radii[row1]
    row2_vdw = vdw_radii[row2]
    vdw_sum = (row1_vdw + row2_vdw) / 10
    vdw_dist = dists <= vdw_sum
    interactions[:, inter_ids['vdw_contacts']] = vdw_dist

    # ---- [HBOND ACCEPTOR] ---------------------------------------------------
    hb_ha_cut = cf.dist_cut_HA
    hb_da_cut = cf.dist_cut_DA
    hb_dha_cut = cf.ang_cut_DHA
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
    angles = calc_angle(xyz2[donors], xyz2[hydros], xyz2[acc])

    count = 0
    for i in range(before_angle.size):
        if before_angle[i]:
            interactions[i, inter_ids['hb_acceptor']] = angles[
                                                            count] > hb_dha_cut
            count += 1

    # ---- [HBOND DONORS] -----------------------------------------------------
    s1_hb_hydros = isin(row1, hb_hydros)
    s2_hb_acc = isin(row2, hb_acc)
    before_angle = hb_dists & s1_hb_hydros & s2_hb_acc

    hydros = row1[before_angle]
    acc = row2[before_angle]
    idx = indices(hb_acc, acc)
    donors = hb_don[idx]
    angles = calc_angle(xyz2[donors], xyz2[hydros], xyz2[acc])

    count = 0
    for i in range(before_angle.size):
        if before_angle[i]:
            interactions[i, inter_ids['hb_donor']] = angles[count] > hb_dha_cut
            count += 1
    return ijf, dists, interactions


max_vdw = iman.get_max_vdw_dist()
to_compute = 'all'
if to_compute == 'all':
    to_compute = interactions
dist_cut = max([y for x in to_compute for y in interactions[x]] + [max_vdw])

s1_indices_raw = iman.sel1_idx
s2_indices_raw = iman.sel2_idx
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
rings = iman.rings

ijf, dists, inters = duplex(xyz, k, inter_ids, s1_indices_raw, s2_indices_raw,
                            dist_cut, anions, cations, hydroph, metal_don,
                            metal_acc, vdw_radii, hb_hydros, hb_donors, hb_acc,
                            rings, to_compute)

# =============================================================================
# %% Start computing interactions
# =============================================================================
# start = time.time()
# import intermap.cutoffs as cf
#
# # ==== close contacts =========================================================
# s1_indices = iman.sel1_idx
# s2_indices = iman.sel2_idx
# dist_cut = cf.close_contacts
# cc_id = inter_identifiers['close_contacts']
# if s1_indices.size and s2_indices.size:
#     cc = close_contacts(xyz, k, cc_id, s1_indices, s2_indices, dist_cut)
#
# # ==== anionic ================================================================
# s1_ani = np.intersect1d(iman.sel1_idx, iman.anions)
# s2_cat = np.intersect1d(iman.sel2_idx, iman.cations)
# anionic_id = inter_identifiers['anionic']
# if s1_ani.size and s2_cat.size:
#     ani = close_contacts(xyz, k, anionic_id, s1_ani, s2_cat, cf.ionic)
#
# # ==== cationic ===============================================================
# s1_cat = np.intersect1d(iman.sel1_idx, iman.cations)
# s2_ani = np.intersect1d(iman.sel2_idx, iman.anions)
# cationic_id = inter_identifiers['cationic']
# if s1_cat.size and s2_ani.size:
#     cat = close_contacts(xyz, k, cationic_id, s1_cat, s2_ani, cf.ionic)
#
# # ==== hydrophobic ============================================================
# s1_hp = np.intersect1d(iman.sel1_idx, iman.hydroph)
# s2_hp = np.intersect1d(iman.sel2_idx, iman.hydroph)
# hydrop_id = inter_identifiers['hydrophobic']
# if s1_hp.size and s2_hp.size:
#     hp = close_contacts(xyz, k, hydrop_id, s1_hp, s2_hp, cf.hydrophobic)
#
# # ==== metal donor ============================================================
# s1_met_donors = np.intersect1d(iman.sel1_idx, iman.metal_don)
# s2_met_acc = np.intersect1d(iman.sel2_idx, iman.metal_acc)
# metd_id = inter_identifiers['metal_donor']
# if s1_met_donors.size and s2_met_acc.size:
#     metd = close_contacts(xyz, k, metd_id, s1_met_donors, s2_met_acc,
#                           cf.metallic)
#
# # ==== metal acceptor =========================================================
# s1_acc = np.intersect1d(iman.sel1_idx, iman.metal_acc)
# s2_met = np.intersect1d(iman.sel2_idx, iman.metal_don)
# meta_id = inter_identifiers['metal_acceptor']
# if s1_acc.size and s2_met.size:
#     meta = close_contacts(xyz, k, meta_id, s1_acc, s2_met, cf.metallic)
#
# # ==== vdw contacts ===========================================================
# max_vdw = iman.get_max_vdw_dist()
# vdw_radii = iman.radii
# vdw_id = inter_identifiers['vdw_contacts']
# if s1_indices.size and s2_indices.size:
#     vdw = close_contacts(xyz, k, vdw_id, s1_indices, s2_indices, max_vdw,
#                          vdw_radii)
#
# # ==== hbonds acceptor ========================================================
# idx = np.isin(iman.hb_H, iman.sel2_idx)
# s2_hb_hydros = iman.hb_H[idx]
# s2_hb_donors = iman.hb_D[idx]
# s1_hb_acc = iman.hb_A[np.isin(iman.hb_A, iman.sel1_idx)]
# hb_acc_id = inter_identifiers['hb_acceptor']
# if s2_hb_donors.size and s1_hb_acc.size:
#     hb_a = dha_contacts(xyz, k, hb_acc_id,
#                         s2_hb_donors, s2_hb_hydros, s1_hb_acc,
#                         cf.HA_cut, cf.DA_cut, cf.DHA_cut)
#
# # ==== hbonds donor ===========================================================
# idx = np.isin(iman.hb_H, iman.sel1_idx)
# s1_hb_hydros = iman.hb_H[idx]
# s1_hb_donors = iman.hb_D[idx]
# s2_hb_acc = iman.hb_A[np.isin(iman.hb_A, iman.sel2_idx)]
# hb_don_id = inter_identifiers['hb_donor']
# if s1_hb_donors.size and s2_hb_acc.size:
#     hb_d = dha_contacts(xyz, k, hb_don_id,
#                         s1_hb_donors, s1_hb_hydros, s2_hb_acc,
#                         cf.HA_cut, cf.DA_cut, cf.DHA_cut)
# print(f"Computing time: {time.time() - start:.2f} s")

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
