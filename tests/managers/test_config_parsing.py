"""
Creator: rglez
Date: 3/22/26
Env:
Desc:
"""
import os
from os.path import join

import pytest

import intermap.managers.config as cfg_mod
from intermap.managers.config import ConfigManager, Param, \
    print_colored_ascii


class TestASCIILogo:
    def test_execution(self, capsys):
        print_colored_ascii()
        assert "Error when reading HTML" not in capsys.readouterr().out

    def test_no_match(self, tmp_path, capsys):
        fake_dir = join(tmp_path, "intermap")
        os.makedirs(fake_dir, exist_ok=True)
        with open(join(fake_dir, "binary_imap.html"), "w") as f:
            f.write("<html>No logo here</html>")

        orig = cfg_mod.proj_dir
        cfg_mod.proj_dir = str(tmp_path)
        try:
            print_colored_ascii()
            assert "\033[38;2;" not in capsys.readouterr().out
        finally:
            cfg_mod.proj_dir = orig

    def test_file_error(self, capsys):
        orig = cfg_mod.proj_dir
        cfg_mod.proj_dir = "/non/existent/path"
        try:
            print_colored_ascii()
            assert "Error when reading HTML" in capsys.readouterr().out
        finally:
            cfg_mod.proj_dir = orig

    def test_empty_lines(self, tmp_path, capsys):
        fake_dir = join(tmp_path, "intermap")
        os.makedirs(fake_dir, exist_ok=True)
        with open(join(fake_dir, "binary_imap.html"), "w") as f:
            f.write('<div style="margin: 20px 0;">Logo\n \nArt</div>')

        orig = cfg_mod.proj_dir
        cfg_mod.proj_dir = str(tmp_path)
        try:
            print_colored_ascii()
            assert "Logo" in capsys.readouterr().out
        finally:
            cfg_mod.proj_dir = orig


class TestBaseParsing:
    def test_param_check_not_implemented(self):
        with pytest.raises(NotImplementedError):
            Param("key", "value").check()

    def test_init_no_path(self, monkeypatch):
        monkeypatch.setattr(cfg_mod, 'detect_config_path',
                            lambda mode: "/fake.cfg")
        with pytest.raises(ValueError, match="No such file"):
            ConfigManager(cfg_path=None)

    def test_corrupt_dtype(self, tmp_path, get_clean_cfg, write_modified_cfg):
        import copy  # <--- Add this
        class BadConfigManager(ConfigManager):
            # Use deepcopy to avoid poisoning the global state!
            allowed_parameters = copy.deepcopy(
                ConfigManager.allowed_parameters)

        BadConfigManager.allowed_parameters['generals']['n_procs'][
            'dtype'] = list
        cfg_dict = get_clean_cfg
        cfg_dict['topo-traj']['topology'] = __file__
        path = write_modified_cfg(tmp_path, cfg_dict, 'generals', 'job_name',
                                  'test')

        with pytest.raises(ValueError, match="dtype is wrong"):
            BadConfigManager(cfg_path=path)


class TestStubbornEdgeCases:
    def test_z_coverage_320_321(self):
        from intermap.managers.config import Config
        c = Config.__new__(Config)
        c.config_path = "dummy"
        c.keyless_sections = []
        c.legal_params = {'sec': {'key': {'dtype': tuple}}}

        class MockConfigObj:
            def sections(self): return ['sec']

            def __getitem__(self, item): return {'key': 'val'}

        c.config_obj = MockConfigObj()
        with pytest.raises(ValueError, match="dtype is wrong"):
            c.check_params()

    def test_z_coverage_445(self):
        from intermap.managers.config import ConfigManager
        cm = ConfigManager.__new__(ConfigManager)
        cm.config_args = {}
        cm.config_obj = {'interactions': {'interactions': 'FakeInteraction'}}

        with pytest.raises(ValueError, match="Invalid interaction"):
            cm.parse_interactions()

    def test_z_coverage_103(self):
        from intermap.managers.config import detect_config_path
        with pytest.raises(ValueError, match="Only modes allowed"):
            detect_config_path(mode='invalid')

    def test_z_coverage_53_to_55(self, tmp_path):
        import intermap.managers.config as cfg
        fake_dir = os.path.join(tmp_path, "intermap")
        os.makedirs(fake_dir, exist_ok=True)
        with open(os.path.join(fake_dir, "binary_imap.html"), "w") as f:
            f.write('<div style="margin: 20px 0;">Logo<br/><br/>Art</div>')

        orig = cfg.proj_dir
        cfg.proj_dir = str(tmp_path)
        try:
            cfg.print_colored_ascii()
        finally:
            cfg.proj_dir = orig
