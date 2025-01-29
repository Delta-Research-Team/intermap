# Created by rglez at 1/19/25
import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd
from numpy import concatenate as concat

import intermap.commons as cmn


@njit(parallel=False, cache=True)
def others(xyz, k, s1_indices, s2_indices, max_vdw, cutoffs_others,
           selected_others):
    """
    Get the containers to find interactions not related to aromatic rings

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_indices (ndarray): Indices of the atoms in the first selection
        s2_indices (ndarray): Indices of the atoms in the second selection
        max_vdw: Maximum van der Waals distance in the current universe
        cutoffs_others (ndarray): Cutoff distances for the interactions
        selected_others (ndarray): Selected nteractions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
        dists (ndarray): Distances between the atoms in the first and second
                         selections
    """

    # Create & query the trees
    dist_cut = max(cutoffs_others[:2].max(), max_vdw)
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], dist_cut)

    # Find the contacts
    n_contacts = sum([len(x) for x in ball_1])
    ijf = np.empty((n_contacts, 3), dtype=np.int32)
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
    dists = cmn.calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    n_types = selected_others.size
    interactions = np.zeros((ijf.shape[0], n_types), dtype=np.bool_)

    return ijf, interactions, dists, row1, row2


@njit(parallel=False, cache=True)
def aro(xyz, k, s1_indices, s2_indices, cations, rings, cutoffs_aro,
        to_compute_aro):
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
        to_compute_aro (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    # Get cations
    s1_cat = s1_indices[cmn.isin(s1_indices, cations)]
    s2_cat = s2_indices[cmn.isin(s2_indices, cations)]

    # Get the aromatic rings
    s1_rings = rings[cmn.isin(rings[:, 0], s1_indices)]
    s2_rings = rings[cmn.isin(rings[:, 0], s2_indices)]

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

    if (s2_indices.size == 0) or (s1_indices.size == 0):
        ijf = np.zeros((0, 3), dtype=np.int32)
        inters = np.zeros((0, len(to_compute_aro)), dtype=np.bool_)
        dists = np.zeros(0, dtype=np.float32)
        row1 = np.zeros(0, dtype=np.int32)
        row2 = np.zeros(0, dtype=np.int32)
        return ijf, inters, dists, row1, row2

    s2_tree = nckd(xyz2[s2_indices])
    ball_1 = s2_tree.query_radius_parallel(xyz2[s1_indices], dist_cut_aro)
    ijf, dists, inters = cmn.get_containers(
        xyz2, k, ext_idx, ball_1, s1_indices, s2_indices, to_compute_aro)
    row1, row2 = ijf[:, 0], ijf[:, 1]
    return ijf, inters, dists, row1, row2
