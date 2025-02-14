# Created by gonzalezroy at 1/23/25
from copy import deepcopy

import pytest

import intermap.config as conf


def save_config(conf_obj, conf_path):
    """
    Save a configuration file
    """
    out_name = conf_path.replace('.cfg', '_mod.cfg')
    with open(out_name, 'w') as configfile:
        conf_obj.write(configfile)
    return out_name


def test_parsing_with_cutoffs_section(conf_obj, conf_path, parameters):
    """Configuration file with a [cutoffs] section is parsed"""
    parsed = conf.ConfigManager(conf_path, parameters)
    assert 'cutoffs' in parsed.config_args


def test_parsing_without_cutoffs_section(conf_obj, conf_path, parameters):
    """
    Configuration file without [cutoffs] section is parsed and that
    section is inferred from internals
    """
    # Save a configuration file without the cutoffs section
    conf_cpy = deepcopy(conf_obj)
    conf_cpy.remove_section('cutoffs')
    out_name = save_config(conf_cpy, conf_path)

    # Parse the modified configuration file
    parsed = conf.ConfigManager(out_name, parameters)
    assert 'cutoffs' in parsed.config_args


def test_invalid_cutoff(conf_obj, conf_path, parameters):
    """
    Invalid cutoff raises an error
    """
    # Save a configuration file with a bad-named cutoff
    conf_cpy = deepcopy(conf_obj)
    conf_cpy.set('cutoffs', 'bad_cutoff', '2')
    out_name = save_config(conf_cpy, conf_path)

    # Parse the modified configuration file
    with pytest.raises(ValueError):
        conf.ConfigManager(out_name, parameters)


def test_parsing_with_missing_cutoffs(conf_obj, conf_path, parameters):
    """
    Configuration file with a [cutoffs] section is parsed
    """
    # Save a configuration file with missing values in [cutoffs] section
    conf_cpy = deepcopy(conf_obj)
    for i, item in enumerate(conf_cpy.items('cutoffs')):
        key, value = item
        if i % 2 == 0:
            conf_cpy.remove_option('cutoffs', key)
    out_name = save_config(conf_cpy, conf_path)

    # Parse the modified configuration file
    parsed = conf.ConfigManager(out_name, parameters)
    assert 'cutoffs' in parsed.config_args

# =============================================================================
#
# =============================================================================
# config = join(conf.proj_dir, 'tests', 'imaps', 'imap1.cfg')
# params = conf.allowed_parameters
# self = TestConfig()
# conf_obj = self.conf_obj(config, params)
