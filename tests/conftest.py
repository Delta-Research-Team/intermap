"""
    Dummy conftest.py for intermap.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    - https://docs.pytest.org/en/stable/fixture.html
    - https://docs.pytest.org/en/stable/writing_plugins.html
"""
import os
from os.path import join

import pytest

from intermap import config as conf


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
    parsed = conf.InterMapConfig(conf_path, parameters)
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
        topo_trajs[case] = (topo, traj)
    return topo_trajs
