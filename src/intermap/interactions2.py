# Created by rglez at 1/19/25
import numpy as np
from numba import njit

from intermap import commons as cmn, containers as cnt


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
    inter_idx = cmn.indices(selected_others, ['VdWContact'])[0]
    vdw_sum = vdw_radii[row1] + vdw_radii[row2]
    return inter_idx, dists <= vdw_sum


@njit(parallel=False, cache=True)
def detect_1d(inter_name, dists, row1, type1, row2, type2, cutoffs_others,
              selected_others):
    """
    Detect the 1D interactions (only one distance involved)
    """

    inter_idx = cmn.indices(selected_others, [inter_name])[0]

    # if type1.size == 0 or type2.size == 0:
    #     return inter_idx, np.zeros(0, dtype=np.bool_)

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
def detect_hbonds(inter_name, row1, type1, row2, type2, dists, xyz, hb_donors,
                  ha_cut, min_ang, max_ang, selected_others):
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
    idx_name = cmn.indices(selected_others, [inter_name])[0]
    r1_t1 = cmn.isin(row1, type1)
    r2_t2 = cmn.isin(row2, type2)
    passing_HA = r1_t1 & r2_t2 & (dists <= ha_cut)
    t1 = row1[passing_HA]
    t2 = row2[passing_HA]
    if not passing_HA.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    # Compute DHA angles
    t3 = np.zeros(passing_HA.size, dtype=np.float32)
    if "HBDonor" in inter_name:
        idx_hydros = cmn.indices(type1, t1)
        D = hb_donors[idx_hydros]
        DHA_angles = cmn.calc_angle(xyz[D], xyz[t1], xyz[t2])
    elif "HBAcceptor" in inter_name:
        idx_hydros = cmn.indices(type2, t2)
        D = hb_donors[idx_hydros]
        DHA_angles = cmn.calc_angle(xyz[D], xyz[t2], xyz[t1])
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")
    t3[passing_HA] = DHA_angles

    # Detect DHA tuples that pass the angle cutoff
    passing_DHA = (t3 >= min_ang) & (t3 <= max_ang)

    if not t3.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)
    return idx_name, passing_DHA


@njit(parallel=False, cache=True)
def detect_2d1a(inter_name, dists, xyz, row1, row2, hb_acc, hb_hydros,
                hb_donors, ha_cut, da_cut, min_ang, max_ang, selected_others):
    """
    Detect the 2D1A interactions (two distances and one angle involved)

    Args:
        inter_name (str): Name of the interaction
        dists (ndarray): Distances between the atoms in the first and second
        xyz (ndarray): Coordinates of the atoms in the system
        row1 (ndarray): Indices of the atoms in the first selection
        row2 (ndarray): Indices of the atoms in the second selection
        hb_acc (ndarray): Indices of the hydrogen bond acceptors
        hb_hydros (ndarray): Indices of the hydrogen bond hydrogens
        hb_donors (ndarray): Indices of the hydrogen bond donors
        ha_cut (float): Cutoff distance for the hydrogen bond acceptor
        da_cut (float): Cutoff distance for the hydrogen bond donor-acceptor
        min_ang (float): Minimum angle for the hydrogen bond
        max_ang (float): Maximum angle for the hydrogen bond
        selected_others (ndarray): Selected interactions to compute

    Returns:
        inter_idx (int): Index of the interaction in the to_compute_others array
        hbonds (ndarray): Boolean array with the atoms that pass the distance
    """
    idx_name = cmn.indices(selected_others, [inter_name])[0]

    if 'Acceptor' in inter_name:
        are_acceptors = cmn.isin(row1, hb_acc)
        are_hydros = cmn.isin(row2, hb_hydros)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        acceptors = row1[passing_ha]
        hydros = row2[passing_ha]

    elif 'Donor' in inter_name:
        are_hydros = cmn.isin(row1, hb_hydros)
        are_acceptors = cmn.isin(row2, hb_acc)
        passing_ha = are_hydros & are_acceptors & (dists <= ha_cut)
        hydros = row1[passing_ha]
        acceptors = row2[passing_ha]
    else:
        raise ValueError(f"Invalid interaction name: {inter_name}")

    if passing_ha.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    ext_idx = cmn.indices(hb_hydros, hydros)
    donors = hb_donors[ext_idx]

    ha_vectors = xyz[acceptors] - xyz[hydros]
    hd_vectors = xyz[donors] - xyz[hydros]
    da_vectors = xyz[donors] - xyz[acceptors]

    da_dists = np.sqrt((da_vectors ** 2).sum(axis=1))
    passing_da = da_dists <= da_cut

    if not passing_da.any():
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    angles = cmn.calc_angles_2v(hd_vectors[passing_da], ha_vectors[passing_da])
    passing_angles = (angles >= min_ang) & (angles <= max_ang)

    if passing_angles.size == 0:
        return idx_name, np.zeros(dists.size, dtype=np.bool_)

    hbonds = np.zeros(dists.size, dtype=np.bool_)
    hbonds[passing_da] = passing_angles
    return idx_name, hbonds


@njit(parallel=False, cache=True)
def fill_others(xyz, k, s1_indices, s2_indices, hydrophobes, anions, cations,
                metal_donors, metal_acceptors, hb_hydros, hb_donors, hb_acc,
                xb_halogens, xb_donors, xb_acc, max_vdw, vdw_radii,
                cutoffs_others, selected_others):
    """
    Fill the not-aromatic interactions
    """

    # Get the containers for the not-aromatic interactions
    ijf, inters, dists, row1, row2 = cnt.others(xyz, k, s1_indices, s2_indices,
                                                max_vdw, cutoffs_others,
                                                selected_others)
    selected = list(selected_others)
    # [Van der Waals]
    if 'VdWContact' in selected:
        vdw_idx, vdw_bit = detect_vdw(dists, row1, row2, vdw_radii,
                                      selected)
        inters[:, vdw_idx] = vdw_bit

    # [Close Contacts]
    if 'CloseContacts' in selected:
        cc_idx, cc = detect_1d('CloseContacts', dists, row1, row1, row2, row2,
                               cutoffs_others, selected)
        inters[:, cc_idx] = cc

    # [Hydrophobic]
    if 'Hydrophobic' in selected:
        hp_idx, hp = detect_1d('Hydrophobic', dists, row1, hydrophobes, row2,
                               hydrophobes, cutoffs_others, selected)
        inters[:, hp_idx] = hp

    # [Cationic]
    if 'Cationic' in selected:
        cat_idx, cat = detect_1d('Cationic', dists, row1, cations, row2,
                                 anions, cutoffs_others, selected)
        inters[:, cat_idx] = cat

    # [Anionic]
    if 'Anionic' in selected:
        an_idx, an = detect_1d('Anionic', dists, row1, anions, row2,
                               cations, cutoffs_others, selected)
        inters[:, an_idx] = an

    # [MetalDonor]
    if 'MetalDonor' in selected:
        md_idx, md = detect_1d('MetalDonor', dists, row1, metal_donors, row2,
                               metal_acceptors, cutoffs_others,
                               selected)
        inters[:, md_idx] = md

    # [MetalAcceptor]
    if 'MetalAcceptor' in selected:
        ma_idx, ma = detect_1d('MetalAcceptor', dists, row1, metal_acceptors,
                               row2, metal_donors, cutoffs_others,
                               selected)
        inters[:, ma_idx] = ma

    # [HBonds]
    if ('HBAcceptor' in selected) or ('HBDonor' in selected):
        idx_hb = cmn.indices(selected, ['HBAcceptor'])[0]
        da_cut_hb, ha_cut_hb, min_ang_hb, max_ang_hb = cutoffs_others[:4,
                                                       idx_hb]
        if 'HBAcceptor' in selected:
            hba_idx, hba = detect_hbonds('HBAcceptor', row1, hb_acc, row2,
                                         hb_hydros, dists, xyz, hb_donors,
                                         ha_cut_hb, min_ang_hb, max_ang_hb,
                                         selected)
            inters[:, hba_idx] = hba

        if 'HBDonor' in selected:
            hbd_idx, hbd = detect_hbonds('HBDonor', row1, hb_hydros, row2,
                                         hb_acc, dists, xyz, hb_donors,
                                         ha_cut_hb, min_ang_hb, max_ang_hb,
                                         selected)
            inters[:, hbd_idx] = hbd

    # [XBonds]
    if ('XBAcceptor' in selected) or ('XBDonor' in selected):
        idx_xb = cmn.indices(selected, ['XBAcceptor'])[0]
        da_cut_xb, ha_cut_xb, min_ang_xb, max_ang_xb = cutoffs_others[:4,
                                                       idx_xb]

        if 'XBAcceptor' in selected:
            xba_idx, xba = detect_hbonds('XBAcceptor', row1, xb_acc, row2,
                                         xb_halogens, dists, xyz, xb_donors,
                                         ha_cut_xb, min_ang_xb, max_ang_xb,
                                         selected)
            inters[:, xba_idx] = xba

        if 'XBDonor' in selected:
            xbd_idx, xbd = detect_hbonds('XBDonor', row1, xb_halogens, row2,
                                         xb_acc, dists, xyz, xb_donors,
                                         ha_cut_xb, min_ang_xb, max_ang_xb,
                                         selected)
            inters[:, xbd_idx] = xbd

    mask = cmn.get_compress_mask(inters)
    return ijf[mask], inters[mask]


# =============================================================================
#
# =============================================================================
from intermap import config as conf
from argparse import Namespace
from intermap.indices import IndexManager
import intermap.cutoffs as cf

conf_path = '/home/rglez/RoyHub/intermap/tests/imaps/imap1.cfg'
# Get the Index Manager
config = conf.InterMapConfig(conf_path, conf.allowed_parameters)
args = Namespace(**config.config_args)
s1 = 'all'
s2 = 'all'
iman = IndexManager(args.topology, args.trajectory, s1, s2, 'all')

# Get information from the Index Manager
u = iman.universe
xyz = u.atoms.positions
k = 0
s1_indices, s2_indices = iman.sel1_idx, iman.sel2_idx
anions, cations = iman.anions, iman.cations
hydrophobes = iman.hydroph
vdw_radii, max_vdw = iman.radii, iman.get_max_vdw_dist()
hb_acc = iman.hb_A
hb_hydros = iman.hb_H
hb_donors = iman.hb_D
xb_acc = iman.xb_A
xb_donors = iman.xb_D
xb_halogens = iman.xb_H
metal_donors = iman.metal_don
metal_acceptors = iman.metal_acc

# Get the interactions and cutoffs
all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
to_compute = all_inters
selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
    cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

ijf, others = fill_others(xyz, k, s1_indices, s2_indices, hydrophobes, anions,
                          cations, metal_donors, metal_acceptors, hb_hydros,
                          hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc,
                          max_vdw, vdw_radii, cutoffs_others, selected_others)
print(others.sum(axis=0))
