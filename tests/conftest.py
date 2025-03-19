"""
    Dummy conftest.py for intermap.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    - https://docs.pytest.org/en/stable/fixture.html
    - https://docs.pytest.org/en/stable/writing_plugins.html
"""
import os
from argparse import Namespace
from os.path import join

import mdtraj as md
import pytest
from interactions import cutoffs as cf
from interactions.indices import IndexManager
from interactions.others import others
from intermap import commons as cmn, config as conf


@pytest.fixture(scope="module")
def conf_path():
    """
    Path to the configuration file for the tests
    """
    return join(conf.proj_dir, 'tests', 'imaps', 'imap1.cfg')


@pytest.fixture(scope="module")
def conf_obj(conf_path, parameters):
    """
    Configuration object for the tests
    """
    parsed = conf.ConfigManager(conf_path, parameters)
    obj = parsed.config_obj
    return obj


@pytest.fixture(scope="module")
def parameters():
    """
    Internally-defined allowed parameters
    """
    return conf.allowed_parameters


@pytest.fixture(scope="module")
def topo_trajs_fix():
    """
    Topology and trajectory files for the tests
    """
    tt_dir = join(conf.proj_dir, 'tests', 'trajs')
    topo_trajs = {}
    for case in os.listdir(tt_dir):
        case_files = os.listdir(join(tt_dir, case))
        topo = None
        traj = None

        for case_file in case_files:
            if case_file.startswith('top'):
                topo = join(tt_dir, case, case_file)
            elif case_file.startswith('traj'):
                traj = join(tt_dir, case, case_file)
            else:
                continue

        assert topo is not None, f"Topology file not found for {case}"
        assert traj is not None, f"Trajectory file not found for {case}"
        if 'trj1' in case:
            topo_trajs[case] = (topo, traj)
    return topo_trajs


@pytest.fixture(scope="module")
def iman(conf_path, parameters):
    """
    Index Manager for the tests
    """
    # Get the Index Manager
    config = conf.ConfigManager(conf_path, parameters)
    args = Namespace(**config.config_args)
    s1 = 'resname LIG'
    s2 = 'protein'
    iman = IndexManager(args.topology, args.trajectory, s1, s2,
                        args.interactions)
    return iman


@pytest.fixture(scope="module")
def others_arguments(conf_path, parameters):
    """
    Arguments for the not-aromatic interactions
    """
    # Get the Index Manager
    config = conf.ConfigManager(conf_path, parameters)
    args = Namespace(**config.config_args)
    s1 = 'resname LIG'
    s2 = 'protein'
    iman = IndexManager(args.topology, args.trajectory, s1, s2,
                        args.interactions)

    # Get information from the Index Manager
    u = iman.universe
    xyz = u.atoms.positions
    k = 0
    s1_indices, s2_indices = iman.sel1_idx, iman.sel2_idx
    anions, cations = iman.anions, iman.cations
    hydrophobes = iman.hydroph
    vdw_radii, max_vdw = iman.radii, iman.get_max_vdw_dist()
    met_donors = iman.metal_don
    met_acc = iman.metal_acc
    hb_hydros = iman.hb_H
    hb_donors = iman.hb_D
    hb_acc = iman.hb_A
    xb_halogens = iman.xb_H
    xb_donors = iman.xb_D
    xb_acc = iman.xb_A

    # Get the interactions and cutoffs
    all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
    to_compute = all_inters
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

    return (
        xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
        met_donors, met_acc, hb_hydros, hb_donors, hb_acc, xb_halogens,
        xb_donors, xb_acc, vdw_radii, max_vdw, cutoffs_others, selected_others,
        args)


@pytest.fixture(scope="module")
def others_arguments_all(conf_path, parameters):
    """
    Arguments for the not-aromatic interactions
    """
    # Get the Index Manager
    config = conf.ConfigManager(conf_path, parameters)
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
    hb_donors = iman.hb_D
    hb_acc = iman.hb_A
    hb_hydros = iman.hb_H
    met_donors = iman.metal_don
    met_acc = iman.metal_acc
    xb_donors = iman.xb_D
    xb_acc = iman.xb_A
    xb_halogens = iman.xb_H

    # Get the interactions and cutoffs
    all_inters, all_cutoffs = cf.get_inters_cutoffs(args.cutoffs)
    to_compute = all_inters
    selected_aro, selected_others, cutoffs_aro, cutoffs_others = \
        cmn.get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs)

    return (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
            vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, met_donors,
            met_acc, xb_donors, xb_acc, xb_halogens, cutoffs_others,
            selected_others, args)


@pytest.fixture(scope="module")
def others_containers(others_arguments):
    """
    Containers for the not-aromatic interactions
    """
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes, met_donors,
     met_acc, hb_hydros, hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc,
     vdw_radii, max_vdw, cutoffs_others, selected_others,
     args) = others_arguments

    ijf, inters = others(xyz, k, s1_indices, s2_indices, hydrophobes, anions,
                         cations, met_donors, met_acc, hb_hydros,
                         hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc,
                         max_vdw, vdw_radii, cutoffs_others, selected_others)
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    return ijf, inters, row1, row2


@pytest.fixture(scope="module")
def others_containers_all(others_arguments_all):
    """
    Containers for the not-aromatic interactions
    """
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     vdw_radii, max_vdw, hb_donors, hb_hydros, hb_acc, met_donors,
     met_acc, xb_donors, xb_acc, xb_halogens, cutoffs_others,
     selected_others, args) = others_arguments_all

    ijf, inters = others(xyz, k, s1_indices, s2_indices, hydrophobes, anions,
                         cations, met_donors, met_acc, hb_hydros,
                         hb_donors, hb_acc, xb_halogens, xb_donors, xb_acc,
                         max_vdw, vdw_radii, cutoffs_others, selected_others)
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    return ijf, inters, row1, row2


@pytest.fixture(scope="module")
def mdtrajectory(others_arguments):
    """
    MDTrajectory object for the tests
    """
    (xyz, k, s1_indices, s2_indices, anions, cations, hydrophobes,
     met_donors, met_acc, hb_hydros, hb_donors, hb_acc, xb_halogens,
     xb_donors, xb_acc, vdw_radii, max_vdw, cutoffs_others, selected_others,
     args) = others_arguments
    trj = md.load(args.trajectory, top=args.topology)[0]
    return trj
