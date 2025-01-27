# Created by rglez at 1/26/25

import MDAnalysis as mda
import pytest
from rgpack.generals import recursive_defaultdict

from intermap import indices as idx


@pytest.fixture(scope="module")
def universes(topo_trajs_fix):
    """
    Universe objects for the tests
    """
    unis = recursive_defaultdict()
    for label, (topo, traj) in topo_trajs_fix.items():
        unis[label]['topo'] = topo
        unis[label]['traj'] = traj
        unis[label]['universe'] = mda.Universe(topo, traj)
    return unis


@pytest.fixture(scope="module")
def iman_universes(universes):
    """
    Universe as parsed by IndexManager
    """
    iman_unis = recursive_defaultdict()
    for label in universes:
        topo = universes[label]['topo']
        traj = universes[label]['traj']
        universe = universes[label]['universe']
        try:
            any_bond = universe.bonds[0]
        except:
            universe = mda.Universe(topo, traj, guess_bonds=True)
            any_bond = universe.bonds[0]
        assert any_bond is not None

        universe.delete_bonds(idx.get_hh_bonds(universe))

        iman_unis[label]['universe'] = universe
        iman_unis[label]['topo'] = topo
        iman_unis[label]['traj'] = traj
    return iman_unis


@pytest.mark.slow
def test_get_uniques(iman_universes):
    """
    """
    for label in iman_universes:
        u = iman_universes[label]['universe']
        unique_mda_res, unique_rdmols, unique_idx = idx.get_uniques(u)
        assert len(unique_mda_res) == len(unique_rdmols) == len(unique_idx)

        if label == 'trj1':
            assert len(unique_mda_res) == 27
            assert sum([len(x) for x in unique_idx]) == 153

        if label == 'trj2':
            assert len(unique_mda_res) == 23
            assert sum([len(x) for x in unique_idx]) == 129

        if label == 'trj3':
            assert len(unique_mda_res) == 37
            assert sum([len(x) for x in unique_idx]) == 209
