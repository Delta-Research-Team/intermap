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
        xyz_chunk, xyzs_aro, xyz_aro_idx, trees_others, trees_aro,
        ijf_shape, inters_shape, s1_idx, s2_idx, anions, cations,
        s1_cat_idx, s2_cat_idx, s1_ani_idx, s2_ani_idx, hydroph, met_don,
        met_acc, vdw_radii, max_vdw, hb_hydr, hb_don, hb_acc, xb_hal, xb_don,
        xb_acc, s1_rings, s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_idx,
        s2_aro_idx, cuts_others, selected_others, cuts_aro, selected_aro,
        overlap, atomic, resconv):
    # Get the number of frames in the chunk
    dim1 = np.int64(xyz_chunk.shape[0])
    dim2 = np.int64(ijf_shape[0])

    # Preallocate the arrays on initialization
    ijf_template = np.empty((dim1, dim2, 3), dtype=np.int32)
    inters_template = np.empty((dim1, dim2, inters_shape[1]), dtype=np.bool_)
    ijf_template.fill(0)
    inters_template.fill(False)

    dist_cut1 = max(cuts_others[:2].max(), max_vdw)
    dist_cut2 = cuts_aro[:2].max()

    # Iterate over the frames in the chunk
    len_aro = len(selected_aro)
    len_others = len(selected_others)
    for i in prange(dim1):
        j = np.int32(i)

        # >>>> Get the coordinates for the frame
        xyz = xyz_chunk[j].copy()
        xyz_aro = xyzs_aro[j].copy()

        # >>>> Get the aromatic interactions
        s1_norm, s2_norm, s1_ctrs, s2_ctrs = aro.get_normals_and_centroids(
            xyz, s1_rings, s2_rings)
        tree_aro = trees_aro[j]
        ball_aro = tree_aro.query_radius(xyz_aro[s1_aro_idx], dist_cut2)

        ijf_aro_tmp, dists_aro = geom.get_containers_run(
            xyz_aro, j, ball_aro, s1_aro_idx, s2_aro_idx)

        ijf_aro, inters_aro = aro.aro(
            xyz_aro, xyz_aro_idx, ijf_aro_tmp, dists_aro, s1_rings_idx,
            s2_rings_idx, s1_cat_idx, s2_cat_idx, s1_ani_idx, s2_ani_idx,
            s1_norm, s2_norm, s1_ctrs, s2_ctrs, cuts_aro, selected_aro)
        num_aros = ijf_aro.shape[0]

        # >>>> Get the non-aromatic interactions
        tree_others = trees_others[j]
        ball_others = tree_others.query_radius(xyz[s1_idx], dist_cut1)

        ijf_others, inters_others = others.others(
            xyz, j, s1_idx, s2_idx, ball_others, hydroph, anions, cations,
            met_don, met_acc, hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc,
            vdw_radii, cuts_others, selected_others, overlap, atomic, resconv)
        num_others = ijf_others.shape[0]

        # >>>> Fill & compress the arrays

        ijf_template[j][:num_aros] = ijf_aro
        ijf_template[j][num_aros:num_aros + num_others] = ijf_others

        inters_template[j][:num_aros, :len_aro] = inters_aro[:num_aros,
                                                  :len_aro]
        inters_template[j][num_aros:num_aros + num_others,
        -len_others:] = inters_others

    ijf_template = ijf_template.reshape(-1, 3)
    inters_template = inters_template.reshape(-1, inters_shape[1])
    mask = geom.get_compress_mask(inters_template)
    return ijf_template[mask], inters_template[mask]


def estimate(
        universe, xyz_aro_idx, chunk_size, s1_idx, s2_idx, cations,
        s1_cat_idx, s2_cat_idx, s1_ani_idx, s2_ani_idx, s1_cat, s2_cat,
        s1_ani, s2_ani, s1_rings, s2_rings, s1_rings_idx, s2_rings_idx,
        s1_aro_idx, s2_aro_idx, cuts_aro, selected_aro, len_aro, anions,
        hydroph, met_don, met_acc, vdw_radii, hb_hydr, hb_don, hb_acc, xb_hal,
        xb_don, xb_acc, cuts_others, selected_others, len_others, max_dist_aro,
        max_dist_others, overlap, atomic, resconv, factor=1.5, n_samples=10):
    # Get the samples coordinates along the trajectory
    n_frames = universe.trajectory.n_frames
    sub = universe.trajectory[::n_frames // n_samples]
    positions = np.asarray([ts.positions.copy() for ts in sub],
                           dtype=np.float32)

    # Detect the number of interactions
    N = len(positions)
    num_detected = np.zeros(N, dtype=np.int32)

    s1_centrs, s2_centrs, xyzs_aro = aro.get_aro_xyzs(
        positions, s1_rings, s2_rings, s1_cat, s2_cat, s1_ani, s2_ani)
    aro_trees = aro.get_trees(xyzs_aro, s2_aro_idx)

    for i in range(N):
        xyz = positions[i]
        xyz_aro = xyzs_aro[i]

        s1_norm, s2_norm, s1_ctrs, s2_ctrs = aro.get_normals_and_centroids(
            xyz, s1_rings, s2_rings)

        tree_aro = aro_trees[i]
        ball_aro = tree_aro.query_radius_parallel(
            xyz_aro[s1_aro_idx], max_dist_aro)

        ijf_aro_tmp, dists, inters = geom.get_containers(
            xyz_aro, i, xyz_aro_idx, ball_aro, s1_aro_idx, s2_aro_idx,
            len(selected_aro), atomic, resconv)

        ijf_aro, inters_aro = aro.aro(
            xyz_aro, xyz_aro_idx, ijf_aro_tmp, dists, s1_rings_idx,
            s2_rings_idx, s1_cat_idx, s2_cat_idx, s1_ani_idx, s2_ani_idx,
            s1_norm, s2_norm, s1_ctrs, s2_ctrs, cuts_aro, selected_aro)

        ball_1 = others.get_ball(xyz, s1_idx, s2_idx, max_dist_others)

        ijf_others, inters_others = others.others(xyz, i, s1_idx, s2_idx,
                                                  ball_1, hydroph, anions,
                                                  cations, met_don, met_acc,
                                                  hb_hydr, hb_don, hb_acc,
                                                  xb_hal, xb_don, xb_acc,
                                                  vdw_radii, cuts_others,
                                                  selected_others, overlap,
                                                  atomic, resconv)
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
