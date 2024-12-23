# Created by rglez at 12/10/24

import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.cutoffs as cf


def parse_cutoff(cutoff_name, args=[]):
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


cutoffs_raw = {
    # Undirected 1D interactions
    'CloseContacts':
        {'distCut1': parse_cutoff('dist_cut_CloseContacts')},
    'VdWContact':
        {'distCut1': 0},
    'Hydrophobic':
        {'distCut1': parse_cutoff('dist_cut_Hydrophobic')},

    # Directed 1D interactions
    'Anionic':
        {'distCut1': parse_cutoff('dist_cut_Ionic')},
    'Cationic':
        {'distCut1': parse_cutoff('dist_cut_Ionic')},
    'MetalDonor':
        {'distCut1': parse_cutoff('dist_cut_Metalic')},
    'MetalAcceptor':
        {'distCut1': parse_cutoff('dist_cut_Metalic')},

    # Directed 2D1A interactions
    'HBAcceptor':
        {'distCut1': parse_cutoff('dist_cut_DA'),
         'distCut2': parse_cutoff('dist_cut_HA')},
    'HBDonor':
        {'distCut1': parse_cutoff('dist_cut_DA'),
         'distCut2': parse_cutoff('dist_cut_HA')},
    'XBAcceptor':
        {'distCut1': parse_cutoff('dist_cut_XA'),
         'distCut2': parse_cutoff('dist_cut_XD')},
    'XBDonor':
        {'distCut1': parse_cutoff('dist_cut_XA'),
         'distCut2': parse_cutoff('dist_cut_XD')},

    # Undirected 2D1A interactions
    'PiStacking':
        {'distCut1': parse_cutoff('dist_cut_PiStacking'),
         'distCut2': parse_cutoff('min_dist_PiStacking'),
         'minAng1': parse_cutoff('min_ang_PiStacking'),
         'maxAng1': parse_cutoff('max_ang_PiStacking')},
    'FaceToFace':
        {'distCut1': parse_cutoff('dist_cut_FaceToFace'),
         'distCut2': parse_cutoff('min_dist_FaceToFace'),
         'minAng1': parse_cutoff('min_ang_FaceToFace'),
         'maxAng1': parse_cutoff('max_ang_FaceToFace')},
    'EdgeToFace':
        {'distCut1': parse_cutoff('dist_cut_EdgeToFace'),
         'distCut2': parse_cutoff('min_dist_EdgeToFace'),
         'minAng1': parse_cutoff('min_ang_EdgeToFace'),
         'maxAng1': parse_cutoff('max_ang_EdgeToFace')},

    # Directed 1D1A interactions
    'PiCation':
        {'distCut1': parse_cutoff('dist_cut_PiCation'),
         'minAng1': parse_cutoff('min_ang_PiCation'),
         'maxAng1': parse_cutoff('max_ang_PiCation')},

    'CationPi':
        {'distCut1': parse_cutoff('dist_cut_PiCation'),
         'minAng1': parse_cutoff('min_ang_PiCation'),
         'maxAng1': parse_cutoff('max_ang_PiCation')}, }


def get_inters_cutoffs(cutoffs_raw):
    """
    Get the interaction names and cutoffs from the raw dictionary

    Args:
        cutoffs_raw (dict): Raw dictionary with the cutoffs.

    Returns:

    """
    inter_names = np.asarray(list(cutoffs_raw.keys()))
    cutoffs_matrix = np.zeros((4, len(inter_names)), dtype=np.float32)

    for i, inter in enumerate(cutoffs_raw):
        if 'distCut1' in cutoffs_raw[inter]:
            cutoffs_matrix[0, i] = cutoffs_raw[inter]['distCut1']

        if 'distCut2' in cutoffs_raw[inter]:
            cutoffs_matrix[1, i] = cutoffs_raw[inter]['distCut2']

        if 'minAng1' in cutoffs_raw[inter]:
            cutoffs_matrix[2, i] = cutoffs_raw[inter]['minAng1']

        if 'maxAng1' in cutoffs_raw[inter]:
            cutoffs_matrix[3, i] = cutoffs_raw[inter]['maxAng1']

    return inter_names, cutoffs_matrix


# =============================================================================
# Helper functions
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
    centr = np.zeros((len(rings), 3), dtype=float)
    for i, ring in enumerate(rings):
        atoms = ring[:ring[-1]]
        centr[i] = xyz[atoms].sum(axis=0) / len(atoms)
    return centr


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


@njit(parallel=False)
def calc_angles_2v(vectors1, vectors2):
    """
    Computes the angles between two arrays of 3D vectors using arctan2.

    Args:
        vectors1 (np.ndarray): Array of shape (N, 3) with the first set of vectors.
        vectors2 (np.ndarray): Array of shape (N, 3) with the second set of vectors.

    Returns:
        np.ndarray: Array of shape (N,) containing the angles in degrees.
    """

    # Input checks and validations
    have_same_shape = vectors1.shape == vectors2.shape
    are_2d_arrays = vectors1.ndim == 2 and vectors2.ndim == 2
    have_3_elements = vectors1.shape[1] == 3
    sms = "Input arrays must be 2D arrays of shape (N, 3) and have the same shape."

    if not (have_same_shape and are_2d_arrays and have_3_elements):
        raise ValueError(sms)
    if vectors1.size == 0:
        return None

    num_vectors = vectors1.shape[0]
    angles_degrees = np.empty(num_vectors)

    for i in range(num_vectors):
        v1 = vectors1[i]
        v2 = vectors2[i]

        dot_product = np.sum(v1 * v2)

        cross_product = np.empty(3)
        cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1]
        cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2]
        cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0]

        cross_product_magnitude = np.sqrt(np.sum(cross_product ** 2))
        angle_radians = np.arctan2(cross_product_magnitude, dot_product)
        angles_degrees[i] = np.degrees(angle_radians)

    return angles_degrees


@njit(parallel=False)
def get_containers(xyz, ext_idx, ball_1, s1_indices, s2_indices, to_compute):
    """
    Get the containers for the interactions, the distances and the indices

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        ext_idx (ndarray): External (real) indices for the atoms in the system
        ball_1 (list): List of lists with the indices of the atoms in the second selection
        s1_indices (ndarray): Indices of the atoms in the first selection
        s2_indices (ndarray): Indices of the atoms in the second selection
        to_compute (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        dists (ndarray): Distances between the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    # Find the contacts
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
    idems = ext_idx[ijf[:, 0]] == ext_idx[ijf[:, 1]]
    if idems.any():
        ijf = ijf[~idems]

    # Compute distances
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    n_types = to_compute.size
    interactions = np.zeros((ijf.shape[0], n_types), dtype=np.bool_)
    return ijf, dists, interactions


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
topo = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
traj = '/media/rglez/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
sel1 = "protein or nucleic or resname Na+"
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

# %% ==========================================================================
# Unified approach
# =============================================================================
import re

inters_all, cutoffs = get_inters_cutoffs(cutoffs_raw)
max_vdw = iman.get_max_vdw_dist()
to_compute = 'all'
if to_compute == 'all':
    to_compute = inters_all

bit_aro = [i for i, x in enumerate(to_compute) if re.search(r'Pi|Face', x)]
bit_others = [i for i in range(cutoffs.shape[1]) if i not in bit_aro]
to_compute_aro = to_compute[bit_aro]
to_compute_others = to_compute[bit_others]

cutoffs_aro = cutoffs[:, bit_aro]
cutoffs_others = cutoffs[:, bit_others]
dist_cut_others = cutoffs_others[[0, 1]].max()

s1_indices_raw = iman.sel1_idx
s2_indices_raw = iman.sel2_idx
vdw_radii = iman.radii
hydroph = iman.hydroph
anions = iman.anions
cations = iman.cations
metal_don = iman.metal_don
metal_acc = iman.metal_acc
hb_hydros = iman.hb_H
hb_donors = iman.hb_D
hb_acc = iman.hb_A
xb_hydros = iman.xb_H
xb_donors = iman.xb_D
xb_acc = iman.xb_A
rings = iman.rings


@njit(parallel=False)
def stackings(inter_name, dists, mindists, s1_normals, s2_normals):
    """
    Helper function to compute the pi-stacking interactions

    """

    # Parse the cutoffs
    idx = indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_dist = cutoffs_aro[1, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Apply restraints
    passing_dist1 = dists <= dist_cut
    passing_dist2 = mindists <= min_dist
    angles = calc_angles_2v(s1_normals, s2_normals)
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
    idx = indices(to_compute_aro, [inter_name])[0]
    dist_cut = cutoffs_aro[0, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Select the pairs
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    if inter_name == 'PiCation':
        s1_is_type = isin(row1, rings_idx)
        s2_is_type = isin(row2, cat_idx)
    elif inter_name == 'CationPi':
        s1_is_type = isin(row1, cat_idx)
        s2_is_type = isin(row2, rings_idx)
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")
    pairs = s1_is_type & s2_is_type

    if pairs.any():
        # Calculate angles between normals and vectors
        row1_pairs = row1[pairs]
        row2_pairs = row2[pairs]
        vector_ctr_cat = xyz2[row1_pairs] - xyz2[row2_pairs]
        normals = rings_normals[indices(rings_idx, row1_pairs)]
        angles = calc_angles_2v(normals, vector_ctr_cat)

        # Apply restraints
        passing_dist = dists[pairs] <= dist_cut
        passing_angles = ((angles >= min_ang) & (angles <= max_ang))
        return idx, pairs, passing_dist & passing_angles
    else:
        return idx, pairs, pairs


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
        ndarray: Container with the aromatic interactions detected
    """
    # =========================================================================
    # STEP I: Find all pair of atoms (or centroids) within the cutoff distance
    # =========================================================================

    # Get cations
    s1_cat = s1_indices_raw[isin(s1_indices_raw, cations)]
    s2_cat = s2_indices_raw[isin(s2_indices_raw, cations)]

    # Get the aromatic rings
    s1_rings = rings[isin(rings[:, 0], s1_indices_raw)]
    s2_rings = rings[isin(rings[:, 0], s2_indices_raw)]

    # Compute the centroids
    s1_centr = compute_centroids(s1_rings, xyz)
    s2_centr = compute_centroids(s2_rings, xyz)

    # Compute the normal vectors
    s1_at1, s1_at3, s1_at5 = s1_rings[:, 0], s1_rings[:, 2], s1_rings[:, 4]
    s2_at1, s2_at3, s2_at5 = s2_rings[:, 0], s2_rings[:, 2], s2_rings[:, 4]
    s1_norm = calc_normal_vector(xyz[s1_at1], xyz[s1_at3], xyz[s1_at5])
    s2_norm = calc_normal_vector(xyz[s2_at1], xyz[s2_at3], xyz[s2_at5])

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

    s2_tree = nckd(xyz2[s2_indices])
    ball_1 = s2_tree.query_radius(xyz2[s1_indices], dist_cut_aro)
    ijf, dists, inters = get_containers(xyz2, ext_idx, ball_1, s1_indices,
                                        s2_indices, to_compute_aro)

    # =========================================================================
    # STEP II: Compute the aromatic interactions
    # =========================================================================
    set_aro = list(to_compute_aro)

    # **** PiCation & CationPi ************************************************
    if 'PiCation' in set_aro:
        idx, pairs, pi_cat = pications('PiCation', ijf, dists, xyz2, s1_norm,
                                       s1_rings_idx, s2_cat_idx,
                                       to_compute_aro)
        inters[pairs, idx] = pi_cat

    if 'CationPi' in set_aro:
        idx, pairs, cat_pi = pications('CationPi', ijf, dists, xyz2, s2_norm,
                                       s2_rings_idx, s1_cat_idx,
                                       to_compute_aro)
        inters[pairs, idx] = cat_pi

    # **** PiStacking, EdgeToFace & FaceToFace ********************************
    find_PiStacking = 'PiStacking' in set_aro
    find_EdgeToFace = 'EdgeToFace' in set_aro
    find_FaceToFace = 'FaceToFace' in set_aro
    if find_PiStacking or find_EdgeToFace or find_FaceToFace:

        # Get the ring pairs
        row1, row2 = ijf[:, 0], ijf[:, 1]
        s1_is_ctr = isin(row1, s1_rings_idx)
        s2_is_ctr = isin(row2, s2_rings_idx)
        pairs = s1_is_ctr & s2_is_ctr
        ring_pairs = ijf[pairs]
        ring_dists = dists[pairs]
        s1_normals = s1_norm[indices(s1_rings_idx, ring_pairs[:, 0])]
        s2_normals = s2_norm[indices(s2_rings_idx, ring_pairs[:, 1])]

        # Compute the minimum distance between the rings
        num_pairs = ring_pairs.shape[0]
        mindists = np.zeros(num_pairs, dtype=np.float32)
        for i in range(num_pairs):
            s1_ring = s1_rings[indices(s1_rings_idx, [ring_pairs[i, 0]])][0]
            s1_ring_idx = s1_ring[:s1_ring[-1]]
            s2_ring = s2_rings[indices(s2_rings_idx, [ring_pairs[i, 1]])][0]
            s2_ring_idx = s2_ring[:s2_ring[-1]]
            mindists[i] = calc_min_dist(xyz[s1_ring_idx], xyz[s2_ring_idx])

        # **** PiStacking *****************************************************
        if find_PiStacking:
            idx, pi_stacking = stackings('PiStacking', ring_dists, mindists,
                                         s1_normals, s2_normals)
            inters[pairs, idx] = pi_stacking

        # **** EdgeToFace *****************************************************
        if find_EdgeToFace:
            idx, etf_stacking = stackings('EdgeToFace', ring_dists, mindists,
                                          s1_normals, s2_normals)
            inters[pairs, idx] = etf_stacking

        # **** FaceToFace *****************************************************
        if find_FaceToFace:
            idx, ftf_stacking = stackings('FaceToFace', ring_dists, mindists,
                                          s1_normals, s2_normals)
            inters[pairs, idx] = ftf_stacking

    return inters


aro_block = aro(xyz, k, s1_indices_raw, s2_indices_raw, cations, rings,cutoffs_aro, to_compute_aro)


# %% ==========================================================================
@njit(parallel=False)
def others(xyz, k, s1_indices_raw, s2_indices_raw, dist_cut, anions,
           cations, hydroph, metal_don, metal_acc, vdw_radii, hb_hydros,
           hb_don, hb_acc, rings, cutoffs, to_compute):
    # =========================================================================
    # STEP I: Find all pair of atoms (and centr) within the cutoff distance
    # =========================================================================

    # Add aromatic centr to xyz
    s1_rings = rings[isin(rings[:, 0], s1_indices_raw)].astype(np.int32)
    s2_rings = rings[isin(rings[:, 0], s2_indices_raw)].astype(np.int32)
    sel1_centr = compute_centr(s1_rings, xyz)
    sel2_centr = compute_centr(s2_rings, xyz)
    xyz2 = np.concatenate((xyz, sel1_centr, sel2_centr), axis=0)

    # Update sel1 & sel2 indices with the aromatic centr
    n0 = xyz.shape[0]
    n1 = n0 + sel1_centr.shape[0]
    n2 = n1 + sel2_centr.shape[0]
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
    inter_ids_selected = [x for x in inter_ids.keys() if x in to_compute]
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
    s1_centr = isin(row1, s1_rings_indices)
    s2_centr = isin(row2, s2_rings_indices)
    close_centr = s1_centr & s2_centr
    ring_pairs = ijf[close_centr]

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
    by_stacking_dist = dists[close_centr] <= stacking_cut
    by_stacking_angles = (angles >= stacking_min) & (angles <= stacking_max)
    by_stacking_mindist = mindists <= stacking_min_dist
    pi_stacking = by_stacking_dist & by_stacking_angles & by_stacking_mindist
    interactions[close_centr, inter_ids['pi_stacking']] = pi_stacking

    # edge-to-face interactions
    by_etf_dist = dists[close_centr] <= etf_cut
    by_etf_angles = (angles >= etf_min) & (angles <= etf_max)
    by_etf_mindist = mindists <= etf_min_dist
    etf = by_etf_dist & by_etf_angles & by_etf_mindist
    interactions[close_centr, inter_ids['edge_to_face']] = etf

    # face-to-face interactions
    by_ftf_dist = dists[close_centr] <= ftf_cut
    by_ftf_angles = (angles >= ftf_min) & (angles <= ftf_max)
    by_ftf_mindist = mindists <= ftf_min_dist
    ftf = by_ftf_dist & by_ftf_angles & by_ftf_mindist
    interactions[close_centr, inter_ids['face_to_face']] = ftf

    # ---- [Pi-Cation] --------------------------------------------------------
    # pi_cat_cut = cf.dist_cut_PiCation
    s2_cat = isin(row2, cations)
    pi_cat = s1_centr & s2_cat
    # pi_cat_dist = dists <= pi_cat_cut
    # compute_angles(s1)
    # interactions[:, inter_ids['pi_cation']] = pi_cat & pi_cat_dist

    # ---- [Cation-Pi] --------------------------------------------------------
    s1_cat = isin(row1, cations)
    # cat_pi = s1_cat & s2_centr
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


ijf, dists, inters = duplex(xyz, k, s1_indices_raw, s2_indices_raw,
                            dist_cut_aro, anions, cations, hydroph, metal_don,
                            metal_acc, vdw_radii, hb_hydros, hb_donors, hb_acc,
                            rings, cutoffs, to_compute)

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
