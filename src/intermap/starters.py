# Created by gonzalezroy at 2/5/25
import logging

import numpy as np
from numba import njit, prange

from intermap.aro import aro
from intermap.others import others

logger = logging.getLogger('InterMapLogger')


def get_estimation(xyz_chunk, n_samples, s1_indices, s2_indices, cations,
                   rings, cutoffs_aro, selected_aro, anions, hydrophobes,
                   metal_donors, metal_acceptors, vdw_radii, max_vdw,
                   hb_hydros, hb_donors, hb_acc, xb_halogens, xb_donors,
                   xb_acc, cutoffs_others, selected_others, factor=2):
    # Preallocate the arrays
    n_frames = xyz_chunk.shape[0]
    samples = xyz_chunk[::n_frames // n_samples]
    num_detected = np.zeros(n_samples, dtype=np.int32)

    # Use parallel loop to fill
    N = samples.shape[0]
    for i in range(N):
        xyz = samples[i]
        ijf_aro, inters_aro = aro(xyz, i, s1_indices, s2_indices, cations,
                                  rings, cutoffs_aro,
                                  selected_aro)

        ijf_others, inters_others = others(xyz, i, s1_indices, s2_indices,
                                           hydrophobes, anions, cations,
                                           metal_donors, metal_acceptors,
                                           hb_hydros, hb_donors, hb_acc,
                                           xb_halogens, xb_donors, xb_acc,
                                           max_vdw, vdw_radii, cutoffs_others,
                                           selected_others)

        num_detected[i] = ijf_aro.shape[0] + ijf_others.shape[0]

    # Estimate the number of contacts
    ijf_vert = np.int32(num_detected.max() * factor)
    if ijf_vert.size == 0:
        ijf_vert = 1500
        logger.warning("No interactions detected in the sampled trajectory."
                       " A default of 1500 cells was automatically allocated ")

    inters_hori = inters_others.shape[1] + inters_aro.shape[1]
    ijf_template = np.zeros((n_frames, ijf_vert, 3), dtype=np.int32)
    inters_template = np.zeros((n_frames, ijf_vert, inters_hori),
                               dtype=np.bool_)
    return ijf_template, inters_template


@njit(parallel=True, cache=True)
def run_parallel(xyz_all, ijf_template, inters_template, len_others, len_aro,
                 s1_indices, s2_indices, anions, cations,
                 hydrophobes, metal_donors, metal_acceptors, vdw_radii,
                 max_vdw, hb_hydros, hb_donors, hb_acc, xb_halogens, xb_donors,
                 xb_acc, rings, cutoffs_others, selected_others, cutoffs_aro,
                 selected_aro):
    num_frames = xyz_all.shape[0]
    limits = np.zeros(num_frames, dtype=np.int32)

    # Use parallel loop to fill
    for i in prange(num_frames):
        xyz = xyz_all[i]

        # Compute the non-aromatic interactions
        ijf_others, inters_others = others(xyz, i, s1_indices, s2_indices,
                                           hydrophobes, anions, cations,
                                           metal_donors, metal_acceptors,
                                           hb_hydros, hb_donors, hb_acc,
                                           xb_halogens, xb_donors, xb_acc,
                                           max_vdw, vdw_radii, cutoffs_others,
                                           selected_others)

        # Compute the aromatic interactions
        ijf_aro, inters_aro = aro(xyz, i, s1_indices, s2_indices, cations,
                                  rings, cutoffs_aro, selected_aro)

        # Fill the templates
        num_others = ijf_others.shape[0]
        num_aro = ijf_aro.shape[0]
        limits[i] = num_others + num_aro

        ijf_template[i, :num_others] = ijf_others
        ijf_template[i, num_others:limits[i]] = ijf_aro
        inters_template[i, :num_others, :len_others] = inters_others
        inters_template[i, num_others:limits[i], len_others:] = inters_aro

    # Compress the chunks to the final size
    num_pairs = limits.sum()
    ijf_final = np.empty((num_pairs, 3), dtype=np.int32)
    inters_final = np.empty((num_pairs, len_others + len_aro), dtype=np.bool_)

    for i in range(num_frames):
        start = limits[:i].sum()
        end = limits[:i + 1].sum()
        ijf_final[start:end] = ijf_template[i, :limits[i]]
        inters_final[start:end] = inters_template[i, :limits[i]]
    return ijf_final, inters_final
