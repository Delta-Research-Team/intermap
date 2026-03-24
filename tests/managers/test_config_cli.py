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

    @pytest.mark.parametrize("mocked_argv", [
        ['intermap', '--help'],  # Triggers help exit
        ['intermap']  # Triggers wrong args exit
    ])
    def test_cli_system_exits(self, monkeypatch, mocked_argv):
        monkeypatch.setattr(sys, 'argv', mocked_argv)
        with pytest.raises(SystemExit):
            detect_config_path(mode='production')

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError, match="Only modes allowed"):
            detect_config_path(mode='unsupported_mode')
