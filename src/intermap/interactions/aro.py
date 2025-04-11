# Created by rglez at 12/10/24
import numpy as np
from numba import njit, prange
from numba.typed import List
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.interactions.geometry as geom


# todo: remove concatenation of arrays and use slicing of preallocated arrays instead

def get_trees(xyzs_aro, s2_aro_indices):
    trees = List()
    for xyz_aro in xyzs_aro:
        s2_tree = nckd(xyz_aro[s2_aro_indices])
        trees.append(s2_tree)
    return trees


@njit(parallel=True, cache=True)
def get_aro_xyzs(xyzs, s1_rings, s2_rings, s1_cat, s2_cat):
    d1, d3 = xyzs.shape[0], xyzs.shape[2]
    d2 = (s1_cat.shape[0] + s2_cat.shape[0] + s1_rings.shape[0] +
          s2_rings.shape[0])

    xyzs_aro = np.empty((d1, d2, d3), dtype=np.float32)
    s1_centrs = np.empty((d1, s1_rings.shape[0], 3), dtype=np.float32)
    s2_centrs = np.empty((d1, s2_rings.shape[0], 3), dtype=np.float32)

    for i in prange(d1):
        xyz = xyzs[i]
        s1_centr = geom.calc_centroids(s1_rings, xyz)
        s2_centr = geom.calc_centroids(s2_rings, xyz)
        xyzs_aro[i] = concat(
            (xyz[s1_cat], xyz[s2_cat], s1_centr, s2_centr), axis=0)
        s1_centrs[i] = s1_centr
        s2_centrs[i] = s2_centr

    return s1_centrs, s2_centrs, xyzs_aro


@njit(parallel=False, cache=True)
def get_normals_and_centroids(xyz, s1_rings, s2_rings):
    s1_at1, s1_at2 = s1_rings[:, 0], s1_rings[:, 1]
    s2_at1, s2_at2 = s2_rings[:, 0], s2_rings[:, 1]
    s1_ctrs = geom.calc_centroids(s1_rings, xyz)
    s2_ctrs = geom.calc_centroids(s2_rings, xyz)
    s1_norm = geom.calc_normal_vector(s1_ctrs, xyz[s1_at1], xyz[s1_at2])
    s2_norm = geom.calc_normal_vector(s2_ctrs, xyz[s2_at1], xyz[s2_at2])
    return s1_norm, s2_norm, s1_ctrs, s2_ctrs


@njit(parallel=False, cache=True)
def get_intersect_point(s1_normal, s1_centroid, s2_normal, s2_centroid):
    # Calculate the intersect direction as the cross product of the two normal vectors
    intersect_direction = np.array([
        s1_normal[1] * s2_normal[2] - s1_normal[2] * s2_normal[1],
        s1_normal[2] * s2_normal[0] - s1_normal[0] * s2_normal[2],
        s1_normal[0] * s2_normal[1] - s1_normal[1] * s2_normal[0],
    ])

    # Set up the system of linear equations to solve
    A = np.stack((s1_normal, s2_normal, intersect_direction))

    # Check if the determinant of A is zero (indicating the vectors are coplanar)
    if np.linalg.det(A) == 0:
        return np.full(3, np.nan, dtype=np.float32)

    # Calculate the offsets for the planes defined by the normals and centroids
    tilted_offset = np.dot(s2_normal, s2_centroid)
    plane_offset = np.dot(s1_normal, s1_centroid)

    # Create the right-hand side of the equation
    d = np.empty(3, dtype=np.float32)
    d[0] = plane_offset
    d[1] = tilted_offset
    d[2] = 0.0

    # Solve the linear equations to find the intersection point on the intersect line
    point = np.linalg.solve(A, d)

    # Find the projection of the centroid onto the intersect line using vector projection
    vec = s1_centroid - point
    intersect_direction_normalized = intersect_direction / np.linalg.norm(
        intersect_direction)
    scalar_proj = np.dot(intersect_direction_normalized, vec)

    return point + (intersect_direction_normalized * scalar_proj)


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
        s1_is_type = geom.isin(row1, s1_rings_idx)
        s2_is_type = geom.isin(row2, s2_cat_idx)
    elif inter_name == 'CationPi':
        s1_is_type = geom.isin(row1, s1_cat_idx)
        s2_is_type = geom.isin(row2, s2_rings_idx)
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
        normals = s1_norm[geom.indices(s1_rings_idx, row1_pairs)]
    else:
        normals = s2_norm[geom.indices(s2_rings_idx, row2_pairs)]

    angles = geom.calc_angles_2v(normals, vector_ctr_cat)

    # Apply restraints
    passing_angles1 = (angles >= min_ang) & (angles <= max_ang)
    passing_angles2 = (angles <= 180 - min_ang) & (angles >= 180 - max_ang)
    passing_angles = passing_angles1 | passing_angles2
    passing = passing_dist & passing_angles
    all_passing = np.zeros(row1.shape[0], dtype=np.bool_)
    all_passing[pairs] = passing
    return idx, all_passing


@njit(parallel=False, cache=True)
def stackings(inter_name, ring_dists, n1n2, n1c1c2, n2c2c1, idists,
              cutoffs_aro, selected_aro):
    """
    Helper function to compute the pi-stacking interactions

    """

    # Parse the cutoffs
    idx = selected_aro.index(inter_name)
    dist_cut = cutoffs_aro[0, idx]
    min_ang1 = cutoffs_aro[2, idx]
    max_ang1 = cutoffs_aro[3, idx]
    min_ang2 = cutoffs_aro[4, idx]
    max_ang2 = cutoffs_aro[5, idx]
    if inter_name == 'FaceToFace':
        iradius = np.inf
    elif inter_name == 'EdgeToFace':
        iradius = 1.5
    elif inter_name == 'PiStacking':
        iradius = np.inf
    else:
        raise Exception(f"Invalid interaction name: {inter_name}")

    # PAssing distances between centroids
    passing_dist = ring_dists <= dist_cut

    # Passing angles between normals
    passing_ang_norms1 = (n1n2 <= max_ang1) & (n1n2 >= min_ang1)
    passing_ang_norms2 = (n1n2 <= 180 - min_ang1) & (n1n2 >= 180 - max_ang1)
    passing_ang_norms = passing_ang_norms1 | passing_ang_norms2

    # Passing angles between normal & centroids
    passing_n1c1c2_1 = (n1c1c2 <= max_ang2) & (n1c1c2 >= min_ang2)
    passing_n1c1c2_2 = (n1c1c2 <= 180 - min_ang2) & (n1c1c2 >= 180 - max_ang2)
    passing_n1c1c2 = passing_n1c1c2_1 | passing_n1c1c2_2

    passing_n2c2c1_1 = (n2c2c1 <= max_ang2) & (n2c2c1 >= min_ang2)
    passing_n2c2c1_2 = (n2c2c1 <= 180 - min_ang2) & (n2c2c1 >= 180 - max_ang2)
    passing_n2c2c1 = passing_n2c2c1_1 | passing_n2c2c1_2

    passing_ncc = passing_n1c1c2 | passing_n2c2c1

    # Passing intersection radius
    passing_radius = idists <= iradius

    stacking = passing_dist & passing_ang_norms & passing_ncc & passing_radius
    return idx, stacking


@njit(parallel=False, cache=True)
def aro(
        xyz_aro, xyz_aro_idx, ijf_aro_tmp, dists_aro, s1_rings_idx,
        s2_rings_idx,
        s1_cat_idx, s2_cat_idx, s1_norm, s2_norm, s1_ctrs, s2_ctrs, cuts_aro,
        selected_aro):
    ijf_aro = ijf_aro_tmp.copy()
    row1 = ijf_aro[:, 0]
    row2 = ijf_aro[:, 1]
    inters_aro = np.zeros((row1.shape[0], len(selected_aro)), dtype=np.bool_)

    if 'None' in selected_aro:
        return ijf_aro, inters_aro

    if 'PiCation' in selected_aro:
        pi_idx, pi_cat = pications(
            'PiCation', xyz_aro, row1, row2, dists_aro, s1_rings_idx,
            s2_rings_idx,
            s1_cat_idx, s2_cat_idx, s1_norm, s2_norm, cuts_aro, selected_aro)
        inters_aro[:, pi_idx] = pi_cat

    if 'CationPi' in selected_aro:
        cat_idx, cat_pi = pications(
            'CationPi', xyz_aro, row1, row2, dists_aro, s1_rings_idx,
            s2_rings_idx,
            s1_cat_idx, s2_cat_idx, s1_norm, s2_norm, cuts_aro, selected_aro)
        inters_aro[:, cat_idx] = cat_pi

    if 'PiStacking' in selected_aro or 'EdgeToFace' in selected_aro or 'FaceToFace' in selected_aro:
        # Get the ring pairs
        s1_is_ctr = geom.isin(row1, s1_rings_idx)
        s2_is_ctr = geom.isin(row2, s2_rings_idx)
        pairs = s1_is_ctr & s2_is_ctr
        ring_pairs = ijf_aro[:, :][pairs]

        s1_idx = geom.indices(s1_rings_idx, ring_pairs[:, 0])
        s2_idx = geom.indices(s2_rings_idx, ring_pairs[:, 1])
        ring_dists = dists_aro[pairs]
        s1_normals, s2_normals = s1_norm[s1_idx], s2_norm[s2_idx]
        s1_centroids, s2_centroids = s1_ctrs[s1_idx], s2_ctrs[s2_idx]

        c1c2 = s1_centroids - s2_centroids
        c2c1 = s2_centroids - s1_centroids
        n1n2 = geom.calc_angles_2v(s1_normals, s2_normals)
        n1c1c2 = geom.calc_angles_2v(c1c2, s1_normals)
        n2c2c1 = geom.calc_angles_2v(c2c1, s2_normals)

        N = s1_normals.shape[0]
        ipoints = np.empty((N, 3), dtype=np.float32)
        for i in range(N):
            s1_normal = s1_normals[i]
            s1_centroid = s1_centroids[i]
            s2_normal = s2_normals[i]
            s2_centroid = s2_centroids[i]
            ipoints[i] = get_intersect_point(s1_normal, s1_centroid,
                                             s2_normal, s2_centroid)

        idists = np.minimum(geom.calc_dist(ipoints, s1_centroids),
                            geom.calc_dist(ipoints, s2_centroids))

        if 'EdgeToFace' in selected_aro:
            idx_etf, etf_stacking = stackings(
                'EdgeToFace', ring_dists, n1n2, n1c1c2, n2c2c1, idists,
                cuts_aro, selected_aro)
            inters_aro[:, idx_etf][pairs] = etf_stacking

        if 'FaceToFace' in selected_aro:
            idx_ftf, ftf_stacking = stackings(
                'FaceToFace', ring_dists, n1n2, n1c1c2, n2c2c1, idists,
                cuts_aro, selected_aro)
            inters_aro[:, idx_ftf][pairs] = ftf_stacking

        if 'PiStacking' in selected_aro:
            idx_pi, pi_stacking = stackings(
                'PiStacking', ring_dists, n1n2, n1c1c2, n2c2c1, idists,
                cuts_aro, selected_aro)
            # inters_aro[:, idx_pi][pairs] = pi_stacking
            inters_aro[:, idx_pi][pairs] = ftf_stacking | etf_stacking

    frame_id = ijf_aro[0][-1]
    mask = geom.get_compress_mask(inters_aro)
    ijf_aro = ijf_aro[mask]
    inters_aro = inters_aro[mask]

    ijf_real_0 = xyz_aro_idx[ijf_aro[:, 0]]
    ijf_real_1 = xyz_aro_idx[ijf_aro[:, 1]]
    ijf_real_2 = np.full(ijf_real_0.shape[0], frame_id, dtype=np.int32)
    ijf_aro = np.stack((ijf_real_0, ijf_real_1, ijf_real_2), axis=1)

    return ijf_aro, inters_aro
