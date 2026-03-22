"""
Creator: rglez
Date: 3/21/26
Env: intermap
Desc: Test fixtures for InterMap.
"""

import configparser
import copy
import os
from os.path import join

import pytest

# Discover paths dynamically
TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
# TESTS_DIR = '/home/rglez/RoyHub/intermap/tests/'
ROOT_DIR = os.path.dirname(TESTS_DIR)

# Point directly to your real config file
CFG_FILE = os.path.join(ROOT_DIR, "examples", "T1", "prot-dna_InterMap.cfg")


@pytest.fixture(scope="session")
def t1_config_dict():
    """
    Reads the real T1 config and returns it as a dictionary.
    """
    if not os.path.exists(CFG_FILE):
        raise FileNotFoundError(
            f"Critical Error: Could not find T1 config at {CFG_FILE}")

    cfg = configparser.ConfigParser()

    cfg.optionxform = str

    cfg.read(CFG_FILE)
    return {section: dict(cfg[section]) for section in cfg.sections()}


@pytest.fixture
def get_clean_cfg(t1_config_dict):
    """
    Provides a fresh, deep-copied configuration dictionary.
    """
    return copy.deepcopy(t1_config_dict)


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
