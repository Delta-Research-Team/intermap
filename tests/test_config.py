# Created by gonzalezroy at 1/23/25
from os.path import join

import pytest

import intermap.config as conf


@pytest.fixture(scope="module")
def config1():
    """
    Path to the configuration file for the tests
    """
    return join(conf.proj_dir, 'tests', 'imaps', 'imap1.cfg')


@pytest.fixture(scope="module")
def allowed_parameters():
    """
    Internally-defined allowed parameters
    """
    return conf.allowed_parameters


class TestConfig:

    def test_overall_parsing(self, config1, allowed_parameters):
        parsed = conf.InterMapConfig(config1, allowed_parameters)
        assert parsed.config_args is not None
