"""
Creator: rglez
Date: 3/22/26
Env:
Desc:
"""
import sys

import pytest

from intermap.managers.config import detect_config_path


class TestConfigCLI:
    def test_production_success(self, monkeypatch):
        monkeypatch.setattr(sys, 'argv', ['intermap', 'my_config.cfg'])
        assert detect_config_path(mode='production') == 'my_config.cfg'

    def test_help_flag_exits(self, monkeypatch):
        monkeypatch.setattr(sys, 'argv', ['intermap', '--help'])
        with pytest.raises(SystemExit):
            detect_config_path(mode='production')

    def test_wrong_args_exits(self, monkeypatch):
        monkeypatch.setattr(sys, 'argv', ['intermap'])
        with pytest.raises(SystemExit):
            detect_config_path(mode='production')

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="Only modes allowed"):
            detect_config_path(mode='unsupported_mode')
