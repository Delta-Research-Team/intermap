# Created by rglez at 1/19/25
"""
This module contains the functions to compute the inter-atomic interactions
"""

import numpy as np
from numba import njit
from numba_kdtree import KDTree as nckd

from intermap.interactions import geometry as aot


def get_ball(xyz, s1_indices, s2_indices, dist_cut):
    s2_tree = nckd(xyz[s2_indices])
    ball_1 = s2_tree.query_radius(xyz[s1_indices], dist_cut)
    return ball_1


@njit(parallel=False, cache=True)
def unswap_frame(ijf):
    """
    Detects swapped pairs in the frame. This occurs when selections overlaps.

    Args:
        ijf: Indices of the atoms in the first and second selections

    Returns:
        uniques: Boolean array to detect the unique pairs
    """
    seen = set()
    uniques = np.zeros(ijf.shape[0], dtype=np.bool_)
    for i in range(ijf.shape[0]):
        a, b, c = ijf[i]
        tupy = (a, b)
        topy = (b, a)
        if (tupy in seen) or (topy in seen):
            continue
        seen.add(tupy)
        seen.add(topy)
        uniques[i] = True

    return uniques


@njit(parallel=False, cache=True)
def containers(xyz, k, s1_idx, s2_idx, ball_1, selected_others,
               overlap, atomic, resconv):
    """
    Get the containers to find interactions not related to aromatic rings

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        s1_idx (ndarray): Indices of the atoms in the first selection
        s2_idx (ndarray): Indices of the atoms in the second selection
        max_vdw: Maximum van der Waals distance in the current universe
        cutoffs_others (ndarray): Cutoff distances for the interactions
        selected_others (ndarray): Selected nteractions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
        dists (ndarray): Distances between the atoms in the first and second
                         selections
    """
    # Find the number of contacts
    n_contacts = sum([len(x) for x in ball_1])

    # Ignore pairs with the same index
    ijf = np.empty((n_contacts, 3), dtype=np.int32)
    counter = 0
    for i, x in enumerate(ball_1):
        X1 = s1_idx[i]
        for j in x:
            X2 = s2_idx[j]
            if X1 == X2:
                continue
            ijf[counter][0] = X1
            ijf[counter][1] = X2
            ijf[counter][2] = k
            counter += 1
    ijf = ijf[:counter]

    # Eliminate the swapped cases if any
    if overlap:
        uniques = unswap_frame(ijf)
        ijf = ijf[uniques]

    # Remove self-residues if not atomic resolution requested
    if not atomic:
        r1 = resconv[ijf[:, 0]]
        r2 = resconv[ijf[:, 1]]
        idems = r1 == r2
        ijf = ijf[~idems]

    # Compute distances
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = aot.calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    N = row1.size
    n_types = len(selected_others)
    interactions = np.zeros((N, n_types), dtype=np.bool_)

    return ijf, interactions, dists, row1, row2


@njit(parallel=False, cache=True)
def detect_vdw(dists, row1, row2, vdw_radii, selected_others):
    """
    Detect the Van der Waals interactions

    Args:
        vdw_radii (ndarray): Van der Waals radii of the atoms
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        dists (ndarray): Distances between the atoms in the first and second
                          selections
        selected_others (ndarray): Selected interactions to compute
    Returns:
        (ndarray): Container with the Van der Waals interactions
    """
    inter_idx = list(selected_others).index('VdWContact')
    vdw_sum = vdw_radii[row1] + vdw_radii[row2]
    return inter_idx, dists <= vdw_sum


@njit(parallel=False, cache=True)
def detect_1d(inter_name, dists, row1, type1, row2, type2, cutoffs_others,
              selected_others):
    """
    Detect the 1D interactions (only one distance involved)
    """
    inter_idx = list(selected_others).index(inter_name)
    dist_cutoff = cutoffs_others[0, inter_idx]
    passing_dists = dists <= dist_cutoff

    if inter_name == 'CloseContact':
        return inter_idx, passing_dists
    else:
        s1_is_type = aot.isin(row1, type1)
        s2_is_type = aot.isin(row2, type2)
        are_type = s1_is_type & s2_is_type
        return inter_idx, passing_dists & are_type


@njit(parallel=False, cache=True)
def detect_hbonds(inter_name, row1, type1, row2, type2, dists, xyz, hb_donors,
                  ha_cut, da_cut, min_ang, max_ang, selected_others):
    """"
    Detect the hydrogen bonds

    Args:
        inter_name (str): Name of the interaction
        row1 (ndarray): Indices of the atoms in the first selection
        type1 (ndarray): Type of the atoms to be found in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        type2 (ndarray): Type of the atoms to be found in the second selection
        dists (ndarray): Distances between the atoms in the first and second
        xyz (ndarray): Coordinates of the atoms in the system
        hb_donors (ndarray): Indices of the hydrogen bond donors
        ha_cut (float): Cutoff distance for the hydrogen-acceptor atoms
        min_ang (float): Minimum angle for the hydrogen bond
        max_ang (float): Maximum angle for the hydrogen bond
        selected_others (ndarray): Selected interactions to compute
    """
    # Detect donors and acceptors under cutoff
    idx_name = list(selected_others).index(inter_name)
    r1_t1 = aot.isin(row1, type1)
    r2_t2 = aot.isin(row2, type2)
    passing_HA = r1_t1 & r2_t2 & (dists <= ha_cut)
    t1 = row1[passing_HA]
    t2 = row2[passing_HA]
    if not passing_HA.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    # Compute DHA angles
    t3 = np.full(passing_HA.size, fill_value=-1, dtype=np.float32)
    if "HBDonor" in inter_name:
        idx_hydros = aot.indices(type1, t1)
        D = hb_donors[idx_hydros]
        DHA_angles = aot.calc_angle_3p(xyz[D], xyz[t1], xyz[t2])
        DA_dists = aot.calc_dist(xyz[D], xyz[t2])
        DHA_angles[DA_dists > da_cut] = -1

    elif "HBAcceptor" in inter_name:
        idx_hydros = aot.indices(type2, t2)
        D = hb_donors[idx_hydros]
        DHA_angles = aot.calc_angle_3p(xyz[D], xyz[t2], xyz[t1])
        DA_dists = aot.calc_dist(xyz[D], xyz[t1])
        DHA_angles[DA_dists > da_cut] = -1
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")
    t3[passing_HA] = DHA_angles

    # Detect DHA tuples that pass the angle cutoff
    passing_DHA1 = (t3 >= min_ang) & (t3 <= max_ang)
    passing_DHA2 = (t3 >= 180 - max_ang) & (t3 <= 180 - min_ang)
    passing_DHA = passing_DHA1 | passing_DHA2

    if not t3.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)
    return idx_name, passing_DHA


@njit(parallel=False, cache=True)
def others(
        xyz, k, s1_idx, s2_idx, ball_others, hydroph, anions, cations, met_don,
        met_acc, hb_hydr, hb_don, hb_acc, xb_hal, xb_don, xb_acc, vdw_radii,
        cuts_others, selected_others, overlap, atomic, resconv):
    """
    Fill the not-aromatic interactions
    """
    if 'None' in selected_others:
        return (np.zeros((0, 3), dtype=np.int32),
                np.zeros((0, 0), dtype=np.bool_))

    # Get the containers for the not-aromatic interactions
    ijf, inters, dists, row1, row2 = \
        containers(xyz, k, s1_idx, s2_idx, ball_others, selected_others,
                   overlap,
                   atomic, resconv)

    # [Van der Waals]
    if 'VdWContact' in selected_others:
        vdw_idx, vdw_bit = detect_vdw(dists, row1, row2, vdw_radii,
                                      selected_others)
        inters[:, vdw_idx] = vdw_bit

    # [Close Contacts]
    if 'CloseContact' in selected_others:
        cc_idx, cc = detect_1d('CloseContact', dists, row1, row1, row2, row2,
                               cuts_others, selected_others)
        inters[:, cc_idx] = cc

    # [Hydrophobic]
    if 'Hydrophobic' in selected_others:
        hp_idx, hp = detect_1d('Hydrophobic', dists, row1, hydroph, row2,
                               hydroph, cuts_others, selected_others)
        inters[:, hp_idx] = hp

    # [Cationic]
    if 'Cationic' in selected_others:
        cat_idx, cat = detect_1d('Cationic', dists, row1, cations, row2,
                                 anions, cuts_others, selected_others)
        inters[:, cat_idx] = cat

    # [Anionic]
    if 'Anionic' in selected_others:
        an_idx, an = detect_1d('Anionic', dists, row1, anions, row2,
                               cations, cuts_others, selected_others)
        inters[:, an_idx] = an

    # [MetalDonor]
    if 'MetalDonor' in selected_others:
        md_idx, md = detect_1d('MetalDonor', dists, row1, met_don, row2,
                               met_acc, cuts_others,
                               selected_others)
        inters[:, md_idx] = md

    # [MetalAcceptor]
    if 'MetalAcceptor' in selected_others:
        ma_idx, ma = detect_1d('MetalAcceptor', dists, row1, met_acc,
                               row2, met_don, cuts_others,
                               selected_others)
        inters[:, ma_idx] = ma

    # [HBonds]
    if ('HBAcceptor' in selected_others) or ('HBDonor' in selected_others):
        idx_hb = selected_others.index('HBAcceptor')
        da_cut_hb, ha_cut_hb, min_ang_hb, max_ang_hb = cuts_others[:4, idx_hb]

        if 'HBAcceptor' in selected_others:
            hba_idx, hba = detect_hbonds('HBAcceptor', row1, hb_acc, row2,
                                         hb_hydr, dists, xyz, hb_don,
                                         ha_cut_hb, da_cut_hb, min_ang_hb,
                                         max_ang_hb, selected_others)
            inters[:, hba_idx] = hba

        if 'HBDonor' in selected_others:
            hbd_idx, hbd = detect_hbonds('HBDonor', row1, hb_hydr, row2,
                                         hb_acc, dists, xyz, hb_don,
                                         ha_cut_hb, da_cut_hb, min_ang_hb,
                                         max_ang_hb, selected_others)
            inters[:, hbd_idx] = hbd

    # [XBonds]
    if ('XBAcceptor' in selected_others) or ('XBDonor' in selected_others):
        idx_xb = selected_others.index('XBAcceptor')
        da_cut_xb, ha_cut_xb, min_ang_xb, max_ang_xb = cuts_others[:4, idx_xb]

        if 'XBAcceptor' in selected_others:
            xba_idx, xba = detect_hbonds('XBAcceptor', row1, xb_acc, row2,
                                         xb_hal, dists, xyz, xb_don,
                                         ha_cut_xb, da_cut_xb, min_ang_xb,
                                         max_ang_xb, selected_others)
            inters[:, xba_idx] = xba

        if 'XBDonor' in selected_others:
            xbd_idx, xbd = detect_hbonds('XBDonor', row1, xb_hal, row2,
                                         xb_acc, dists, xyz, xb_don,
                                         ha_cut_xb, da_cut_xb, min_ang_xb,
                                         max_ang_xb, selected_others)
            inters[:, xbd_idx] = xbd

    mask = aot.get_compress_mask(inters)
    return ijf[mask], inters[mask]
