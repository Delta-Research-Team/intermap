# Created by rglez at 1/28/25

import mdtraj as md
import numpy as np

import intermap.commons as cmn
import intermap.interactions2 as int2


# from tests.conftest import mdtrajectory, others_containers

class TestGetIndexInter:
    def test_get_index_inter_found(self):
        """Returns correct index when interaction is found"""
        inter_name = 'Hydrophobic'
        selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
        result = cmn.indices([inter_name], selected_others)
        position = np.where(result != -1)[0][0]
        assert position == 1

    def test_get_index_inter_not_found(self):
        """Returns -1 when interaction is not found"""
        inter_name = 'Anionic'
        selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
        result = cmn.indices(inter_name, selected_others)
        position = result[result != -1]
        assert position.size == 0

    def test_get_index_inter_empty_list(self):
        """Returns -1 when selected_others list is empty"""
        inter_name = 'Hydrophobic'
        selected_others = np.array([])
        result = cmn.indices(inter_name, selected_others)
        position = result[result != -1]
        assert position.size == 0

    def test_get_index_inter_multiple_occurrences(self):
        """Returns the first index when interaction appears multiple times"""
        inter_name = 'Hydrophobic'
        selected_others = np.array(
            ['VdWContact', 'Hydrophobic', 'Hydrophobic'])
        result = cmn.indices([inter_name], selected_others)
        position = np.where(result != -1)[0][0]
        assert position == 1


def test_distances(mdtrajectory, others_containers):
    ijf, _, dists, _, _ = others_containers
    dists_mdtraj = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    are_equals = np.isclose(dists, dists_mdtraj)
    assert are_equals.size == are_equals.sum()


def test_vdw(others_arguments_all, others_containers_all, mdtrajectory):
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, cutoffs_others,
     selected_others, args) = others_arguments_all

    ijf, _, dists, row1, row2 = others_containers_all
    n, vdw = int2.detect_vdw(dists, row1, row2, vdw_radii, selected_others)

    dists_mdtraj = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    r1 = vdw_radii[ijf[:, 0]]
    r2 = vdw_radii[ijf[:, 1]]
    contacts = dists_mdtraj <= (r1 + r2)
    assert np.all(contacts == vdw)


def test_close_contacts(others_containers_all, others_arguments_all,
                        mdtrajectory):
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, cutoffs_others,
     selected_others, args) = others_arguments_all

    ijf, _, dists, row1, row2 = others_containers_all
    n, cc = int2.detect_1d('CloseContacts', dists, row1, row1, row2, row2,
                           cutoffs_others, selected_others)
    dist_cutoff = cutoffs_others[0, n]

    # Test all imap cc are under the cutoff
    imap_dists = md.compute_distances(mdtrajectory, ijf[cc][:, :2]) * 10
    assert imap_dists.max() <= dist_cutoff

    # Test all the distances
    md_dists = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    assert np.all(np.isclose(md_dists, dists))


def test_1d(others_containers_all, others_arguments_all, mdtrajectory):
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, cutoffs_others,
     selected_others, args) = others_arguments_all

    ijf, inters, dists, row1, row2 = others_containers_all

    inter_idx = {'Hydrophobic': {'r1': hydrophobes, 'r2': hydrophobes},
                 'Anionic': {'r1': anions, 'r2': cations},
                 'Cationic': {'r1': cations, 'r2': anions}}

    for inter_name in inter_idx:
        type1 = inter_idx[inter_name]['r1']
        type2 = inter_idx[inter_name]['r2']
        n, inter_bit = int2.detect_1d(inter_name, dists, row1, type1, row2,
                                      type2, cutoffs_others, selected_others)

        if inter_bit.sum() == 0:
            continue

        dist_cutoff = cutoffs_others[0, n]

        md_dists = md.compute_distances(mdtrajectory,
                                        ijf[:, :2][inter_bit]) * 10
        assert np.all(md_dists <= dist_cutoff)


def test_hbonds(others_containers_all, others_arguments_all, mdtrajectory):
    # Arrange data
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, cutoffs_others,
     selected_others, args) = others_arguments_all
    ijf, _, dists, row1, row2 = others_containers_all

    # Get the cutoffs for the hydrogen bonds
    ha_cut = 2.5
    min_ang = 120
    max_ang = 180

    # Detect hbonds in mdtraj (using Hubbard's method)
    hubards = md.baker_hubbard(mdtrajectory)[:, [1, 2]]
    unique_ha_hubard = np.unique(hubards, axis=0)
    unique_ha_hubard_tuples = tuple(map(tuple, unique_ha_hubard))

    # Detect hbonds in intermap
    hb1 = int2.detect_hbonds.py_func('HBDonor', row1, hb_hydros, row2, hb_acc,
                                     dists, xyz, hb_donors, ha_cut, min_ang,
                                     max_ang, selected_others)

    unique_ha_iman = np.unique(ijf[hb1[1]][:, :2], axis=0)
    unique_ha_iman_tuples = tuple(map(tuple, unique_ha_iman))

    # Test all Hubbard's hbonds are in intermap's hbonds
    for tuple_hubard in unique_ha_hubard_tuples:
        assert tuple_hubard in unique_ha_iman_tuples, f"Tuple {tuple_hubard} not found"

    # Test with mdtraj that all distances of imap's hbonds are under the cutoff
    distances = md.compute_distances(mdtrajectory,
                                     unique_ha_iman_tuples) * 10
    assert np.all(distances <= ha_cut), f"Distances not met"

    # Test with mdtraj that all angles of imap's hbonds are under the cutoff
    H = ijf[hb1[1]][:, 0]
    D = hb_donors[cmn.indices(hb_hydros, H)]
    A = ijf[hb1[1]][:, 1]
    DHA = np.concatenate((D, H, A), axis=0).reshape(3, -1).T
    angles = md.compute_angles(mdtrajectory, DHA) * 180 / np.pi
    assert np.all((angles >= min_ang) & (angles <= max_ang))
