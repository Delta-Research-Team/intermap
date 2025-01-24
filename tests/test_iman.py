# Created by gonzalezroy at 1/24/25
from os.path import join

import pytest

import intermap.config as conf
from intermap import indices as idx


@pytest.fixture(scope="module")
def iman():
    """
    Path to a topology file for the tests
    """
    config_path = join(conf.proj_dir, 'tests', 'imaps', 'imap1.cfg')
    parsed = conf.InterMapConfig(config_path, conf.allowed_parameters)
    topo = parsed.config_args['topology']
    traj = parsed.config_args['trajectory']
    inters = parsed.config_args['interactions']
    sel1 = parsed.config_args['selection_1']
    sel2 = parsed.config_args['selection_2']
    return idx.IndexManager(topo, traj, sel1, sel2, inters)


class TestIman:

    def test_hh(self, iman):
        """
        """
        universe = iman.universe
        assert not idx.any_hh_bonds(universe)

    def test_get_hh_bonds(self, iman):
        """
        """
        universe = iman.universe
        with pytest.raises(StopIteration):
            next(idx.get_hh_bonds(universe))

    def test_get_uniques(self, iman):
        """
        """
        universe = iman.universe
        unique_mda_res, unique_rdmols, unique_idx = idx.get_uniques(universe)
        assert len(unique_mda_res) == len(unique_rdmols) == len(
            unique_idx) == 27
        assert sum([len(x) for x in unique_idx]) == 153
