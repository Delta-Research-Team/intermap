# Created by rglez at 1/19/25
import numpy as np
from numba import njit

import intermap.commons as cmn
import intermap.containers as cnt


# todo: remove bonded from the interactions

@njit(parallel=False, cache=True)
def get_index_inter(inter_name, selected_others):
    """
    Get the index of an interaction in the selected_others list
    """
    where = np.where(selected_others == inter_name)[0]
    if where.size == 0:
        return -1
    return where[0]


@njit(parallel=False, cache=True)
def detect_vdw(inter_name, dists, row1, row2, vdw_radii, selected_others):
    """
    Detect the Van der Waals interactions

    Args:
        inter_name (str): Name of the interaction
        vdw_radii (ndarray): Van der Waals radii of the atoms
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        dists (ndarray): Distances between the atoms in the first and second
                          selections
        selected_others (ndarray): Selected interactions to compute
    Returns:
        (ndarray): Container with the Van der Waals interactions
    """
    inter_idx = cmn.indices(selected_others, [inter_name])[0]
    vdw_sum = vdw_radii[row1] + vdw_radii[row2]
    return inter_idx, dists <= vdw_sum


@njit(parallel=False, cache=True)
def detect_1d(inter_name, dists, row1, type1, row2, type2, cutoffs_others,
              selected_others):
    """
    Detect the 1D interactions (only one distance involved)
    """

    inter_idx = cmn.indices(selected_others, [inter_name])[0]

    if type1.size == 0 or type2.size == 0:
        return inter_idx, np.zeros(0, dtype=np.bool_)

    dist_cutoff = cutoffs_others[0, inter_idx]
    passing_dists = dists <= dist_cutoff

    if inter_name == 'CloseContacts':
        return inter_idx, passing_dists
    else:
        s1_is_type = cmn.isin(row1, type1)
        s2_is_type = cmn.isin(row2, type2)
        are_type = s1_is_type & s2_is_type
        return inter_idx, passing_dists & are_type


@njit(parallel=False, cache=True)
def fill_others(xyz, k, s1_indices, s2_indices, hydrophobes, anions, cations,
                max_vdw, vdw_radii, cutoffs_others, selected_others):
    """
    Fill the not-aromatic interactions
    """
    selected = list(selected_others)
    ijf, inters, dists, row1, row2 = cnt.others(xyz, k, s1_indices, s2_indices,
                                                max_vdw, cutoffs_others,
                                                selected_others)

    if 'VdWContact' in selected:
        vdw_index, vdw_contacts = detect_vdw('VdWContact', dists, row1, row2,
                                             vdw_radii, selected_others)
        inters[:, vdw_index] = vdw_contacts

    if 'CloseContacts' in selected:
        cc_index, close_contacts = detect_1d('CloseContacts', dists, row1,
                                             row1, row2, row2, cutoffs_others,
                                             selected_others)
        inters[:, cc_index] = close_contacts

    if 'Hydrophobic' in selected:
        hp_index, hp_contacts = detect_1d('Hydrophobic', dists, row1,
                                          hydrophobes, row2, hydrophobes,
                                          cutoffs_others, selected_others)
        inters[:, hp_index] = hp_contacts

    exist_salt_bridges = (anions.size > 0) and (cations.size > 0)
    find_cationic = ('Cationic' in selected) and exist_salt_bridges
    find_anionic = ('Anionic' in selected) and exist_salt_bridges

    if find_cationic:
        cat_idx, cationic = detect_1d('Cationic', dists, row1, cations, row2,
                                      anions, cutoffs_others, selected_others)
        inters[:, cat_idx] = cationic

    if find_anionic:
        cat_idx, anionic = detect_1d('Anionic', dists, row1, anions, row2,
                                     cations, cutoffs_others, selected_others)
        inters[:, cat_idx] = anionic
