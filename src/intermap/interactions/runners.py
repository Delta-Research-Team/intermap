# Created by gonzalezroy at 2/5/25
import logging
from sys import getsizeof

import numpy as np
from numba import njit, prange

import intermap.interactions.geometry as geom
from intermap.interactions import aro, others

logger = logging.getLogger('InterMapLogger')


@njit(parallel=True, cache=True)
def runpar(
        xyz_chunk, xyzs_aro, xyz_aro_real_idx, trees_others, trees_aro,
        ijf_shape, inters_shape, s1_indices, s2_indices, anions, cations,
        s1_cat_idx, s2_cat_idx, hydrophobes, metal_donors, metal_acceptors,
        vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc, xb_halogens,
        xb_donors, xb_acc, s1_rings, s2_rings, s1_rings_idx, s2_rings_idx,
        s1_aro_indices, s2_aro_indices, cutoffs_others, selected_others,
        cutoffs_aro, selected_aro, overlap):
    # Get the number of frames in the chunk
    dim1 = np.int64(xyz_chunk.shape[0])
    dim2 = np.int64(ijf_shape[0])

    # Preallocate the arrays on initialization
    ijf_template = np.empty((dim1, dim2, 3), dtype=np.int32)
    inters_template = np.empty((dim1, dim2, inters_shape[1]), dtype=np.bool_)
    ijf_template.fill(0)
    inters_template.fill(False)

    markers = np.zeros(dim1, dtype=np.int32)
    dist_cut1 = max(cutoffs_others[:2].max(), max_vdw)
    dist_cut2 = cutoffs_aro[:2].max()
    n_aros = np.zeros(dim1, dtype=np.int32)
    n_others = np.zeros(dim1, dtype=np.int32)

    # Iterate over the frames in the chunk
    for i in prange(dim1):
        # Get the coordinates for the frame
        xyz = xyz_chunk[i]
        xyz_aro = xyzs_aro[i]
        ijf_view = ijf_template[i].copy()
        inters_view = inters_template[i].copy()

        # >>>> Get the aromatic interactions
        tree_aro = trees_aro[i]
        ball_aro = tree_aro.query_radius(xyz_aro[s1_aro_indices], dist_cut2)

        s1_norm, s2_norm, s1_ctrs, s2_ctrs = aro.get_normals_and_centroids(
            xyz, s1_rings, s2_rings)
        ijf_view, dists = geom.get_containers_run(
            xyz_aro, i, xyz_aro_real_idx, ball_aro, s1_aro_indices,
            s2_aro_indices, ijf_view)
        n_pairs_aro = dists.shape[0]
        ijf_view, inters_view, n_aro = aro.aro(xyz_aro, xyz_aro_real_idx,
                                               n_pairs_aro, ijf_view, dists,
                                               inters_view, s1_rings_idx,
                                               s2_rings_idx, s1_cat_idx,
                                               s2_cat_idx, s1_norm, s2_norm,
                                               s1_ctrs, s2_ctrs, cutoffs_aro,
                                               selected_aro)
        n_aros[i] = n_aro

        # >>>> Get the non-aromatic interactions
        tree_others = trees_others[i]
        ball_others = tree_others.query_radius(xyz[s1_indices], dist_cut1)
        ijf_others, inters_others = others.others(
            xyz, i, s1_indices, s2_indices, ball_others, hydrophobes, anions,
            cations, metal_donors, metal_acceptors, hb_hydros, hb_donors,
            hb_acc, xb_halogens, xb_donors, xb_acc, vdw_radii, cutoffs_others,
            selected_others, overlap)
        n_others[i] = ijf_others.shape[0]

        M = inters_view.shape[1] - inters_others.shape[1]
        ijf_view[n_aros[i]:n_aros[i] + n_others[i]] = ijf_others
        inters_view[n_aros[i]:n_aros[i] + n_others[i], M:] = inters_others
        markers[i] = n_aros[i] + n_others[i]

        ijf_template[i] = ijf_view
        inters_template[i] = inters_view

    ijf_reshaped = ijf_template.reshape(-1, 3)
    inters_reshaped = inters_template.reshape(-1, inters_shape[1])
    mask = geom.get_compress_mask(inters_reshaped)
    return ijf_reshaped[mask], inters_reshaped[mask]


def estimate(
        universe, xyz_aro_idx, chunk_size, s1_idx, s2_idx, cations,
        s1_cat_idx, s2_cat_idx, s1_cat, s2_cat, s1_rings, s2_rings,
        s1_rings_idx, s2_rings_idx, s1_aro_idx, s2_aro_idx, cuts_aro,
        selected_aro, len_aro, anions, hydroph, met_don, met_acc, vdw_radii,
        hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc, cuts_others,
        selected_others, len_others, max_dist_aro, max_dist_others, overlap,
        atomic, resconv, factor=1.5, n_samples=10):
    # Get the samples coordinates along the trajectory

    n_frames = universe.trajectory.n_frames
    sub = universe.trajectory[::n_frames // n_samples]
    positions = np.asarray([ts.positions.copy() for ts in sub],
                           dtype=np.float32)

    # Detect the number of interactions
    N = len(positions)
    num_detected = np.zeros(N, dtype=np.int32)

    s1_centrs, s2_centrs, xyzs_aro = aro.get_aro_xyzs(
        positions, s1_rings, s2_rings, s1_cat, s2_cat)
    aro_trees = aro.get_trees(xyzs_aro, s2_aro_idx)
    for i in range(N):
        xyz = positions[i]
        xyz_aro = xyzs_aro[i]

        s1_norm, s2_norm, s1_ctrs, s2_ctrs = aro.get_normals_and_centroids(
            xyz, s1_rings, s2_rings)

        tree_aro = aro_trees[i]
        ball_aro = tree_aro.query_radius_parallel(xyz_aro[s1_aro_idx],
                                                  max_dist_aro)

        ijf, dists, inters = geom.get_containers(xyz_aro, i, xyz_aro_idx,
                                                 ball_aro, s1_aro_idx,
                                                 s2_aro_idx, len(selected_aro),
                                                 atomic, resconv)

        n_pairs_aro = dists.shape[0]

        ijf_aro, inters_aro, n_aro = aro.aro(xyz_aro, xyz_aro_idx, n_pairs_aro,
                                             ijf, dists, inters, s1_rings_idx,
                                             s2_rings_idx, s1_cat_idx,
                                             s2_cat_idx, s1_norm, s2_norm,
                                             s1_ctrs, s2_ctrs, cuts_aro,
                                             selected_aro)

        ball_1 = others.get_ball(xyz, s1_idx, s2_idx, max_dist_others)

        ijf_others, inters_others = others.others(
            xyz, i, s1_idx, s2_idx, ball_1, hydroph, anions, cations, met_don,
            met_acc, hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc,
            vdw_radii, cuts_others, selected_others, overlap)
        num_detected[i] = ijf_aro.shape[0] + ijf_others.shape[0]

    # Estimate the number of contacts
    ijf_vert = np.int32(num_detected.max() * factor)
    if ijf_vert < 1000:
        ijf_vert = 1500
    inters_hori = len_others + len_aro

    # Get the sizes for preallocation
    ijf_dummy = np.zeros((ijf_vert, np.int32(3)))
    ijf_size = round(getsizeof(ijf_dummy) * chunk_size / 2 ** 20)
    inters_dummy = np.zeros((ijf_vert, inters_hori), dtype=np.bool_)
    inters_size = round(getsizeof(inters_dummy) * chunk_size / 2 ** 20)

    n_coords = positions[0].shape[0] * chunk_size + xyzs_aro[0].shape[
        0] * chunk_size
    coords_dummy = np.zeros((n_coords, 3), dtype=np.float32)
    coords_size = round(getsizeof(coords_dummy) / 2 ** 20)

    round((getsizeof(positions[0]) * chunk_size + getsizeof(
        xyzs_aro[0]) * chunk_size) / 2 ** 20)

    logger.info(
        f"Estimated base memory allocation for coordinates "
        f"({coords_size} MB) and interactions ({ijf_size + inters_size} MB) arrays.")

    ijf_shape = ijf_dummy.shape
    inters_shape = inters_dummy.shape
    return ijf_shape, inters_shape
