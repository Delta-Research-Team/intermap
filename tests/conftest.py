"""
Creator: rglez
Date: 3/21/26
Env: intermap
Desc: Test fixtures for InterMap.
"""

import configparser
import copy
import os
from os.path import abspath, dirname, exists, join
from types import SimpleNamespace

import numpy as np
import pytest
import scipy.spatial

import intermap.commons as cmn
import intermap.interactions.aro as aro_module
import intermap.interactions.others as others_module
from intermap.managers.indices import get_periodic_table_info, IndexManager

# =============================================================================
# PATHS & CONFIG
# =============================================================================
TESTS_DIR = dirname(abspath(__file__))
ROOT_DIR = dirname(TESTS_DIR)
CFG_FILE = join(ROOT_DIR, "examples", "T1", "prot-dna_InterMap.cfg")


# =============================================================================
# FIXTURES
# =============================================================================


@pytest.fixture(autouse=True)
def fix_numba_recursion(monkeypatch):
    """
    Swaps Numba-KDTree for SciPy-KDTree to avoid RecursionError
    while maintaining real scientific data flow.
    """

    class ScipyKDTree:
        def __init__(self, data, *args, **kwargs):
            self.tree = scipy.spatial.cKDTree(data)

        def query_radius(self, x, r):
            return [np.array(i, dtype=np.int64) for i in
                    self.tree.query_ball_point(x, r)]

        def query_radius_parallel(self, x, r):
            return self.query_radius(x, r)

    # Patch the direct local references where KDTree was imported as 'nckd'
    monkeypatch.setattr(cmn, "nckd", ScipyKDTree)
    monkeypatch.setattr(aro_module, "nckd", ScipyKDTree)
    monkeypatch.setattr(others_module, "nckd", ScipyKDTree)


@pytest.fixture(scope="session")
def t1_config_dict():
    """
    Reads the real T1 config and returns it as a dictionary.
    """
    if not exists(CFG_FILE):
        raise FileNotFoundError(
            f"Critical Error: Could not find T1 config at {CFG_FILE}")

    cfg = configparser.ConfigParser()

    cfg.optionxform = str

    cfg.read(CFG_FILE)
    return {section: dict(cfg[section]) for section in cfg.sections()}


@pytest.fixture(scope="session")
def t2_system_args():
    """Arguments for the T2 system (OPC water + Lone Pairs)."""
    base_dir = abspath(join(dirname(__file__), 'data', 'T2'))
    return SimpleNamespace(
        topology=join(base_dir, 'generic.pdb'),
        trajectory=join(base_dir, 'generic.dcd'),
        selection_1="protein",
        selection_2="water",
        interactions="all",
        start=0,
        last=2,
        stride=1,
        resolution='atom',
        output_dir=join(base_dir, 'temp_test_out_t2'),
        annotations=False,
        cutoffs={},
        chunk_size=1,
        job_name="T2_test",
        min_prevalence=0.0  # <--- ADDED
    )


@pytest.fixture(scope="session")
def example_system_args():
    """Provides a Namespace mimicking ConfigManager for the T1 example data."""
    base_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', 'examples', 'T1'))
    return SimpleNamespace(
        topology=os.path.join(base_dir,
                              'hmr_cionize_ions_solvent_sc150mM_WRAPPED_RENUM.pdb'),
        trajectory=os.path.join(base_dir, 'traj_sample.dcd'),
        selection_1="protein",
        selection_2="nucleic",
        interactions="all",
        start=0,
        last=2,
        stride=1,
        resolution='atom',
        output_dir=os.path.abspath(
            os.path.join(os.path.dirname(__file__), 'temp_test_out')),
        annotations=False,
        cutoffs={},
        chunk_size=1,
        job_name="T1_test",
        min_prevalence=0.0  # <--- ADDED
    )


@pytest.fixture
def get_clean_cfg(t1_config_dict):
    """
    Provides a fresh config with safe dummy paths already set.
    """
    cfg = copy.deepcopy(t1_config_dict)
    # Set safe defaults so validation doesn't trip on file checks
    cfg['topo-traj']['topology'] = __file__
    cfg['topo-traj']['trajectory'] = "dummy"
    return cfg


@pytest.fixture
def write_modified_cfg():
    """
    Returns a function to write a modified config file preserving case.
    """

    def _write(tmp_path, base_dict, section, key, new_value):
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        for s, params in base_dict.items():
            cfg.add_section(s)
            for k, v in params.items():
                cfg[s][k] = str(v)
        if new_value is None:
            if key in cfg[section]:
                del cfg[section][key]
        else:
            if not cfg.has_section(section):
                cfg.add_section(section)
            cfg[section][key] = str(new_value)
        cfg_path = join(tmp_path, f"mod_{section}_{key}.cfg")
        with open(cfg_path, "w") as f:
            cfg.write(f)
        return str(cfg_path)

    return _write


@pytest.fixture(scope="session")
def iman(example_system_args, tmp_path_factory):
    """
    Loads IndexManager once per session.
    """
    # Create a unique temp directory for this session's output
    temp_out = str(tmp_path_factory.mktemp("shared_iman_out"))
    example_system_args.output_dir = temp_out
    return IndexManager(example_system_args)


@pytest.fixture(scope="session")
def iman_t2(t2_system_args, tmp_path_factory):  # Removed 'self'
    """
    Loads IndexManager for T2 system once per session.
    """

    temp_out = str(tmp_path_factory.mktemp("shared_t2_out"))
    # Ensure we use a fresh copy of args
    args = copy.deepcopy(t2_system_args)
    args.output_dir = temp_out
    return IndexManager(args)


@pytest.fixture(scope="session")
def pt_info():
    """
    Provides the periodic table dictionary used for element guessing.
    """
    return get_periodic_table_info()
