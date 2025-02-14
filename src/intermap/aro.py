# Created by rglez at 12/10/24
import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.geometry as aot


# todo: remove concatenation of arrays and use slicing of preallocated arrays instead


@njit(parallel=False, cache=True)
def containers_aro(xyz, k, s1_indices, s2_indices, cations, rings, cutoffs_aro,
                   selected_aro):
    """
    Compute the aromatic interactions

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices (ndarray): Indices of the atoms in the first selection
        s2_indices (ndarray): Indices of the atoms in the second selection
        cations (ndarray): Indices of the cations
        rings (ndarray): Indices of the aromatic rings
        cutoffs_aro (ndarray): Cutoff distances for the aromatic interactions
        selected_aro (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    # Get cations
    s1_cat = s1_indices[aot.isin(s1_indices, cations)]
    s2_cat = s2_indices[aot.isin(s2_indices, cations)]

    # Get the aromatic rings
    s1_rings = rings[aot.isin(rings[:, 0], s1_indices)]
    s2_rings = rings[aot.isin(rings[:, 0], s2_indices)]

    # Compute the centroids
    s1_centr = aot.calc_centroids(s1_rings, xyz)
    s2_centr = aot.calc_centroids(s2_rings, xyz)

    # Compute the normal vectors
    s1_at1, s1_at3, s1_at5 = s1_rings[:, 0], s1_rings[:, 2], s1_rings[:, 4]
    s2_at1, s2_at3, s2_at5 = s2_rings[:, 0], s2_rings[:, 2], s2_rings[:, 4]
    s1_norm = aot.calc_normal_vector(xyz[s1_at1], xyz[s1_at3], xyz[s1_at5])
    s2_norm = aot.calc_normal_vector(xyz[s2_at1], xyz[s2_at3], xyz[s2_at5])

    # Create a new xyz array with the cations & centroids only
    xyz_aro = concat((xyz[s1_cat], xyz[s2_cat], s1_centr, s2_centr), axis=0)
    xyz_aro_real_idx = concat(
        (s1_cat, s2_cat, s1_rings[:, 0], s2_rings[:, 0])).astype(np.int32)

    # Internal indexing for xyz2 coordinates
    n0 = s1_cat.size + s2_cat.size
    n1 = n0 + s1_centr.shape[0]
    n2 = n1 + s2_centr.shape[0]

    s1_cat_idx = np.arange(0, s1_cat.size, dtype=np.int32)
    s2_cat_idx = np.arange(s1_cat.size, n0, dtype=np.int32)
    s1_rings_idx = np.arange(n0, n1, dtype=np.int32)
    s2_rings_idx = np.arange(n1, n2, dtype=np.int32)

    # Create & query the trees
    dist_cut = cutoffs_aro[:2].max()
    s1_aro_indices = concat((s1_cat_idx, s1_rings_idx)).astype(np.int32)
    s2_aro_indices = concat((s2_cat_idx, s2_rings_idx)).astype(np.int32)
    s2_tree = nckd(xyz_aro[s2_aro_indices])
    ball_1 = s2_tree.query_radius(xyz_aro[s1_aro_indices], dist_cut)

    # Get containers
    ijf, dists, interactions = aot.get_containers(
        xyz_aro, k, xyz_aro_real_idx, ball_1, s1_aro_indices, s2_aro_indices,
        len(selected_aro))

    row1, row2 = ijf[:, 0], ijf[:, 1]

    return (xyz_aro, xyz_aro_real_idx, s1_cat_idx, s2_cat_idx, s1_rings_idx,
            s2_rings_idx, ijf, interactions, dists, row1, row2, s1_norm,
            s2_norm, s1_rings, s2_rings)


@njit(parallel=False, cache=True)
def pications(inter_name, xyz_aro, row1, row2, dists, s1_rings_idx,
              s2_rings_idx, s1_cat_idx, s2_cat_idx, s1_norm, s2_norm,
              cutoffs_aro, selected_aro):
    """
    Helper function to compute the pi-cation // cation-pi interactions

    """
    # Parse the cutoffs
    idx = selected_aro.index(inter_name)
    dist_cut = cutoffs_aro[0, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Select the pairs
    if inter_name == 'PiCation':
        s1_is_type = aot.isin(row1, s1_rings_idx)
        s2_is_type = aot.isin(row2, s2_cat_idx)
    elif inter_name == 'CationPi':
        s1_is_type = aot.isin(row1, s1_cat_idx)
        s2_is_type = aot.isin(row2, s2_rings_idx)
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")

    pairs = s1_is_type & s2_is_type
    passing_dist = dists[pairs] <= dist_cut
    if not passing_dist.any():
        return idx, np.zeros(row1.shape[0], dtype=np.bool_)

    # Calculate angles between normals and vectors
    row1_pairs = row1[pairs]
    row2_pairs = row2[pairs]

    vector_ctr_cat = xyz_aro[row1_pairs] - xyz_aro[row2_pairs]

    if inter_name == 'PiCation':
        normals = s1_norm[aot.indices(s1_rings_idx, row1_pairs)]
    else:
        normals = s2_norm[aot.indices(s2_rings_idx, row2_pairs)]

    angles = aot.calc_angles_2v(normals, vector_ctr_cat)

    # Apply restraints
    passing_angles = ((angles >= min_ang) & (angles <= max_ang))
    passing = passing_dist & passing_angles
    all_passing = np.zeros(row1.shape[0], dtype=np.bool_)
    all_passing[pairs] = passing
    return idx, all_passing


@njit(parallel=False, cache=True)
def stackings(inter_name, ring_dists, mindists, s1_normals, s2_normals,
              cutoffs_aro, selected_aro):
    """
    Helper function to compute the pi-stacking interactions

    """

    # Parse the cutoffs
    idx = selected_aro.index(inter_name)
    dist_cut = cutoffs_aro[0, idx]
    min_dist = cutoffs_aro[1, idx]
    min_ang = cutoffs_aro[2, idx]
    max_ang = cutoffs_aro[3, idx]

    # Apply restraints
    passing_dist1 = ring_dists <= dist_cut
    passing_dist2 = (mindists <= min_dist) & (mindists > 0)
    angles = aot.calc_angles_2v(s1_normals, s2_normals)
    passing_angles = ((angles >= min_ang) & (angles <= max_ang))
    stacking = passing_dist1 & passing_dist2 & passing_angles
    return idx, stacking


# @njit("Tuple((i4[:, :], b1[:, :]))(f4[:, :], i4, i4[:], i4[:], i4[:], i4[:, :], i8[:], i8[:])", parallel=False, cache=True)
@njit(parallel=False, cache=True)
def aro(xyz, k, s1_indices, s2_indices, cations, rings, cutoffs_aro,
        selected_aro):
    """
    Compute the aromatic interactions

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        cations (ndarray): Indices of the cations
        rings (ndarray): Indices of the aromatic rings
        cutoffs_aro (ndarray): Cutoff distances for the aromatic interactions
        selected_aro (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """
    if 'None' in selected_aro:
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, 0), dtype=np.bool_)
        return ijf, inters
    # Get containers
    (xyz_aro, xyz_aro_real_idx, s1_cat_idx, s2_cat_idx, s1_rings_idx,
     s2_rings_idx, ijf, inters, dists, row1, row2, s1_norm,
     s2_norm, s1_rings, s2_rings) = containers_aro(xyz, k, s1_indices,
                                                   s2_indices, cations, rings,
                                                   cutoffs_aro, selected_aro)

    if 'PiCation' in selected_aro:
        pi_idx, pi_cat = pications('PiCation', xyz_aro, row1, row2, dists,
                                   s1_rings_idx, s2_rings_idx, s1_cat_idx,
                                   s2_cat_idx, s1_norm, s2_norm, cutoffs_aro,
                                   selected_aro)
        inters[:, pi_idx] = pi_cat

    if 'CationPi' in selected_aro:
        cat_idx, cat_pi = pications('CationPi', xyz_aro, row1, row2, dists,
                                    s1_rings_idx, s2_rings_idx, s1_cat_idx,
                                    s2_cat_idx, s1_norm, s2_norm, cutoffs_aro,
                                    selected_aro)
        inters[:, cat_idx] = cat_pi

    if 'PiStacking' in selected_aro or 'EdgeToFace' in selected_aro or 'FaceToFace' in selected_aro:
        # Get the ring pairs
        s1_is_ctr = aot.isin(row1, s1_rings_idx)
        s2_is_ctr = aot.isin(row2, s2_rings_idx)
        pairs = s1_is_ctr & s2_is_ctr

        ring_pairs = ijf[pairs]
        ring_dists = dists[pairs]
        s1_normals = s1_norm[aot.indices(s1_rings_idx, ring_pairs[:, 0])]
        s2_normals = s2_norm[aot.indices(s2_rings_idx, ring_pairs[:, 1])]

        # Compute the minimum distance between the rings
        num_pairs = ring_pairs.shape[0]
        mindists = np.zeros(num_pairs, dtype=np.float32)
        for i in range(num_pairs):
            r1 = (s1_rings[:, 0] == xyz_aro_real_idx[
                ring_pairs[i, 0]]).argmax()
            ring1 = s1_rings[r1][:s1_rings[r1][-1]]
            r2 = (s2_rings[:, 0] == xyz_aro_real_idx[
                ring_pairs[i, 1]]).argmax()
            ring2 = s2_rings[r2][:s2_rings[r2][-1]]
            mindists[i] = aot.calc_min_dist(xyz[ring1], xyz[ring2])

        if 'PiStacking' in selected_aro:
            idx, pi_stacking = stackings('PiStacking', ring_dists, mindists,
                                         s1_normals, s2_normals, cutoffs_aro,
                                         selected_aro)

            inters[pairs, idx] = pi_stacking

        if 'EdgeToFace' in selected_aro:
            idx, etf_stacking = stackings('EdgeToFace', ring_dists, mindists,
                                          s1_normals, s2_normals, cutoffs_aro,
                                          selected_aro)
            inters[pairs, idx] = etf_stacking

        if 'FaceToFace' in selected_aro:
            idx, ftf_stacking = stackings('FaceToFace', ring_dists, mindists,
                                          s1_normals, s2_normals, cutoffs_aro,
                                          selected_aro)
            inters[pairs, idx] = ftf_stacking

    mask = aot.get_compress_mask(inters)
    ijf_mask = ijf[mask]
    inters_mask = inters[mask]
    ijf_real_row1 = xyz_aro_real_idx[ijf_mask[:, 0]]
    ijf_real_row2 = xyz_aro_real_idx[ijf_mask[:, 1]]
    ijf_real_row3 = ijf_mask[:, 2]
    ijf_real = np.stack((ijf_real_row1, ijf_real_row2, ijf_real_row3), axis=1)
    return ijf_real, inters_mask
