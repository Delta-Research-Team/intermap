# Created by rglez at 1/28/25

import mdtraj as md
import numpy as np

import intermap.commons as cmn
import intermap.interactions2 as int2
# from tests.conftest import mdtrajectory, others_containers


def test_get_index_inter_found():
    """Returns correct index when interaction is found"""
    inter_name = 'Hydrophobic'
    selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
    result = cmn.indices([inter_name], selected_others)
    position = np.where(result != -1)[0][0]
    assert position == 1


def test_get_index_inter_not_found():
    """Returns -1 when interaction is not found"""
    inter_name = 'Anionic'
    selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
    result = cmn.indices(inter_name, selected_others)
    position = result[result != -1]
    assert position.size == 0


def test_get_index_inter_empty_list():
    """Returns -1 when selected_others list is empty"""
    inter_name = 'Hydrophobic'
    selected_others = np.array([])
    result = cmn.indices(inter_name, selected_others)
    position = result[result != -1]
    assert position.size == 0


def test_get_index_inter_multiple_occurrences():
    """Returns the first index when interaction appears multiple times"""
    inter_name = 'Hydrophobic'
    selected_others = np.array(['VdWContact', 'Hydrophobic', 'Hydrophobic'])
    result = cmn.indices([inter_name], selected_others)
    position = np.where(result != -1)[0][0]
    assert position == 1


def test_distances(mdtrajectory, others_containers):
    ijf, _, dists, _, _ = others_containers
    dists_mdtraj = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    are_equals = np.isclose(dists, dists_mdtraj)
    assert are_equals.size == are_equals.sum()


def test_vdw(others_arguments, others_containers, mdtrajectory):
    _, _, _, _, _, _, _, vdw_radii, _, _, selected_others, _ = others_arguments

    ijf, _, dists, row1, row2 = others_containers
    n, vdw = int2.detect_vdw(dists, row1, row2, vdw_radii, selected_others)

    dists_mdtraj = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    r1 = vdw_radii[ijf[:, 0]]
    r2 = vdw_radii[ijf[:, 1]]
    contacts = dists_mdtraj <= (r1 + r2)
    assert np.all(contacts == vdw)


def test_close_contacts(others_containers, others_arguments, mdtrajectory):
    (_, _, s1_indices, s2_indices, _, _, _, _, _, cutoffs_others,
     selected_others, _) = others_arguments

    ijf, _, dists, row1, row2 = others_containers
    n, cc = int2.detect_1d('CloseContacts', dists, row1, row1, row2, row2,
                           cutoffs_others, selected_others)
    dist_cutoff = cutoffs_others[0, n]

    # Test the total number of contacts
    product = cmn.product_uniques(s1_indices, s2_indices)
    md_dists_all = md.compute_distances(mdtrajectory, product) * 10
    n_cc = (md_dists_all <= dist_cutoff).sum()
    assert n_cc == cc.sum()

    # Test the individual contacts
    md_dists = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
    md_cc = md_dists <= dist_cutoff
    assert np.all(md_cc == cc)


def test_1d(others_containers, others_arguments, mdtrajectory):
    (_, _, _, _, anions, cations, hydrophobes, _, _, cutoffs_others,
     selected_others, _) = others_arguments

    ijf, inters, dists, row1, row2 = others_containers

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

        s1_idx = row1[cmn.isin(row1, type1)]
        s2_idx = row2[cmn.isin(row2, type2)]

        product = cmn.product_uniques(s1_idx, s2_idx)
        md_dists_all = md.compute_distances(mdtrajectory, product) * 10
        n_inters = (md_dists_all <= dist_cutoff).sum()
        assert n_inters == inter_bit.sum()

        md_dists = md.compute_distances(mdtrajectory, ijf[:, :2]) * 10
        assert np.all(md_dists <= dist_cutoff)
