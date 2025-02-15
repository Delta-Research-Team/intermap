# Created by gonzalezroy at 2/5/25
import logging

import numpy as np
from numba import njit, prange

import intermap.commons as cmn
import intermap.geometry as geom
from intermap import aro, others

logger = logging.getLogger('InterMapLogger')


@njit(parallel=True, cache=True)
def runpar(xyz_chunk, xyzs_aro, xyz_aro_real_idx, trees_chunk, aro_balls,
           ijf_shape, inters_shape, len_others, len_aro, s1_indices,
           s2_indices, anions, cations, s1_cat_idx, s2_cat_idx, hydrophobes,
           metal_donors, metal_acceptors, vdw_radii, max_vdw, hb_hydros,
           hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc, s1_rings,
           s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_indices,
           s2_aro_indices, cutoffs_others, selected_others, cutoffs_aro,
           selected_aro):
    # Get the number of frames in the chunk
    dim1 = np.int64(xyz_chunk.shape[0])
    dim2 = np.int64(ijf_shape[0])

    # Preallocate the arrays
    ijf_template = np.zeros((dim1, dim2, 3), dtype=np.int32)
    inters_template = np.zeros((dim1, dim2, inters_shape[1]), dtype=np.bool_)
    markers = np.zeros(dim1, dtype=np.int32)
    dist_cut1 = max(cutoffs_others[:2].max(), max_vdw)

    for i in prange(dim1):
        # Compute the interactions
        xyz = xyz_chunk[i]
        xyz_aro = xyzs_aro[i]
        s1_norm, s2_norm = aro.get_aro_normals(xyz, s1_rings, s2_rings)
        ball_aro = aro_balls[i]

        ijf, dists, inters = geom.get_containers(
            xyz_aro, i, xyz_aro_real_idx, ball_aro, s1_aro_indices,
            s2_aro_indices, len(selected_aro))

        ijf_aro, inters_aro = aro.aro(
            xyz, xyz_aro, xyz_aro_real_idx, ijf, dists, inters, s1_rings,
            s1_rings_idx, s2_rings, s2_rings_idx, s1_cat_idx, s2_cat_idx,
            s1_norm, s2_norm, cutoffs_aro, selected_aro)

        tree_others = trees_chunk[i]
        ball_others = cmn.get_ball(xyz, s1_indices, tree_others, dist_cut1)
        ijf_others, inters_others = others.others(
            xyz, i, s1_indices, s2_indices, ball_others, hydrophobes, anions,
            cations, metal_donors, metal_acceptors, hb_hydros, hb_donors,
            hb_acc, xb_halogens, xb_donors, xb_acc, vdw_radii, cutoffs_others,
            selected_others)

        # Fill the templates
        num_others = len(ijf_others)
        num_aro = len(ijf_aro)
        markers[i] = num_others + num_aro

        ijf_template[i][:num_others] = ijf_others
        ijf_template[i][num_others:markers[i]] = ijf_aro
        inters_template[i][:num_others, :len_others] = inters_others
        inters_template[i][num_others:markers[i], len_others:] = inters_aro

    # Compress the chunks to the final size
    num_pairs = markers.sum()
    ijf_final = np.empty((num_pairs, 3), dtype=np.int32)
    inters_final = np.empty((num_pairs, len_others + len_aro), dtype=np.bool_)

    for i in prange(dim1):
        start = markers[:i].sum()
        end = markers[:i + 1].sum()
        ijf_final[start:end] = ijf_template[i, :markers[i]]
        inters_final[start:end] = inters_template[i, :markers[i]]
    return ijf_final, inters_final


def estimate(
        positions, xyz_aro_real_idx, chunk_size, s1_indices, s2_indices,
        cations, s1_cat_idx, s2_cat_idx, s1_cat, s2_cat, s1_rings,
        s2_rings, s1_rings_idx, s2_rings_idx, s1_aro_indices, s2_aro_indices,
        cutoffs_aro, selected_aro, len_aro, anions, hydrophobes, metal_donors,
        metal_acceptors, vdw_radii, max_vdw, hb_hydros, hb_donors, hb_acc,
        xb_halogens, xb_donors, xb_acc, cutoffs_others, selected_others,
        len_others, dist_cut_aro, factor=1.5):

    # Detect the number of interactions
    N = len(positions)
    others_cut = max(cutoffs_others[:2].max(), max_vdw) if len_others else 0
    num_detected = np.zeros(N, dtype=np.int32)

    s1_centrs, s2_centrs, xyzs_aro = aro.get_aro_xyzs(
        positions, s1_rings, s2_rings, s1_cat, s2_cat)
    aro_balls = aro.get_balls(xyzs_aro, s1_aro_indices, s2_aro_indices,
                              dist_cut_aro)

    for i in prange(N):
        xyz = positions[i]
        xyz_aro = xyzs_aro[i]
        s1_norm, s2_norm = aro.get_aro_normals(xyz, s1_rings, s2_rings)
        ball_aro = aro_balls[i]
        ijf, dists, inters = geom.get_containers(
            xyz_aro, i, xyz_aro_real_idx, ball_aro, s1_aro_indices,
            s2_aro_indices, len(selected_aro))
        ijf_aro, inters_aro = aro.aro(xyz, xyz_aro, xyz_aro_real_idx, ijf,
                                      dists, inters,
                                      s1_rings, s1_rings_idx, s2_rings,
                                      s2_rings_idx, s1_cat_idx, s2_cat_idx,
                                      s1_norm, s2_norm, cutoffs_aro,
                                      selected_aro)
        ball_1 = others.get_ball(xyz, s1_indices, s2_indices, others_cut)
        ijf_others, inters_others = others.others(xyz, i, s1_indices,
                                                  s2_indices, ball_1,
                                                  hydrophobes, anions,
                                                  cations,
                                                  metal_donors,
                                                  metal_acceptors, hb_hydros,
                                                  hb_donors, hb_acc,
                                                  xb_halogens, xb_donors,
                                                  xb_acc, vdw_radii,
                                                  cutoffs_others,
                                                  selected_others)

        num_detected[i] = ijf_aro.shape[0] + ijf_others.shape[0]

    # Estimate the number of contacts
    ijf_vert = np.int32(num_detected.max() * factor)
    if ijf_vert == 0:
        ijf_vert = 1500
    inters_hori = len_others + len_aro

    # Get the shapes for preallocation
    ijf_shape = (ijf_vert, np.int32(3))
    inters_shape = (ijf_vert, inters_hori)

    # Get the numbers in MB
    v_size = ijf_shape[0]
    h_size = ijf_shape[1] + inters_shape[1]
    ijf_mb = chunk_size * v_size * 3 * 32 / 2 ** 20
    inters_mb = chunk_size * v_size * h_size * 8 / 2 ** 20
    mb1 = np.float32(round(ijf_mb + inters_mb, 0))
    mb2 = np.float32(round(chunk_size * (
        np.union1d(s1_indices, s2_indices).size) * 3 * 4 / 2 ** 20, 0))
    return ijf_shape, inters_shape, mb1, mb2, v_size, h_size
