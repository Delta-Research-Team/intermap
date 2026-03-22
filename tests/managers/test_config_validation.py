"""
Creator: rglez
Date: 3/22/26
Env:
Desc:
"""
import os
from os.path import dirname, exists, join, normpath

import pytest

from intermap.managers.config import ConfigManager


class TestDirectoryAndPaths:
    def test_output_dir_normalization(self, tmp_path, get_clean_cfg,
                                      write_modified_cfg):
        cfg = get_clean_cfg
        cfg['topo-traj']['topology'] = __file__
        cfg['topo-traj']['trajectory'] = "dummy"
        path = write_modified_cfg(tmp_path, cfg, 'generals', 'output_dir',
                                  "results")

        cm = ConfigManager(cfg_path=path)
        assert cm.config_args['output_dir'] == normpath(
            join(dirname(path), "results"))

    def test_directory_creation(self, tmp_path, get_clean_cfg,
                                write_modified_cfg):
        cfg = get_clean_cfg
        out_dir = join(tmp_path, "new_out")
        cfg['topo-traj']['topology'] = __file__
        cfg['topo-traj']['trajectory'] = "dummy"
        cfg['generals']['output_dir'] = out_dir

        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name',
                                  'test_job')
        ConfigManager(cfg_path=path)
        assert exists(out_dir)
        assert exists(join(out_dir, "test_job_InterMap.cfg"))

    def test_directory_already_exists(self, tmp_path, get_clean_cfg,
                                      write_modified_cfg):
        cfg = get_clean_cfg
        out_dir = join(tmp_path, "already_here")
        os.makedirs(out_dir, exist_ok=True)
        cfg['topo-traj']['topology'] = __file__
        cfg['generals']['output_dir'] = out_dir

        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name',
                                  'test')
        with pytest.raises(FileExistsError):
            ConfigManager(cfg_path=path)


class TestParameterConstraints:
    def test_missing_section(self, tmp_path, get_clean_cfg):
        import configparser
        cfg_obj = configparser.ConfigParser()
        for s, params in get_clean_cfg.items():
            if s != 'generals': cfg_obj[s] = params

        path = join(tmp_path, "missing_sec.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)
        with pytest.raises(ValueError, match="Missing sections"):
            ConfigManager(cfg_path=path)

    @pytest.mark.parametrize("val, match", [("0", "out of"), ("-1", "out of"),
                                            ("abc", "int")])
    def test_n_procs_validation(self, tmp_path, get_clean_cfg,
                                write_modified_cfg, val, match):
        path = write_modified_cfg(tmp_path, get_clean_cfg, 'generals',
                                  'n_procs', val)
        with pytest.raises(ValueError, match=match):
            ConfigManager(cfg_path=path)

    def test_bad_cutoff(self, tmp_path, get_clean_cfg, write_modified_cfg):
        cfg = get_clean_cfg
        cfg['topo-traj']['topology'] = __file__
        cfg['topo-traj']['trajectory'] = "dummy"
        cfg['generals']['output_dir'] = join(tmp_path, "err_cut")
        cfg['cutoffs']['fake_cutoff'] = '5.0'

        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name', 'c')
        with pytest.raises(ValueError, match="is not a valid cutoff name"):
            ConfigManager(cfg_path=path)

    def test_bad_interaction(self, tmp_path, get_clean_cfg,
                             write_modified_cfg):
        cfg = get_clean_cfg
        cfg['topo-traj']['topology'] = __file__
        cfg['topo-traj']['trajectory'] = "dummy"
        cfg['generals']['output_dir'] = join(tmp_path, "err_int")

        path = write_modified_cfg(tmp_path, cfg, 'interactions',
                                  'interactions', 'MagicForce')
        with pytest.raises(ValueError, match="Invalid interaction"):
            ConfigManager(cfg_path=path)


class TestDroppedCoverage:
    def test_missing_cutoffs_section(self, tmp_path, get_clean_cfg):
        import configparser
        cfg_dict = get_clean_cfg
        cfg_dict['topo-traj']['topology'] = __file__
        cfg_dict['topo-traj']['trajectory'] = "dummy"

        cfg_obj = configparser.ConfigParser()
        cfg_obj.optionxform = str
        for s, p in cfg_dict.items():
            if s != 'cutoffs':
                cfg_obj[s] = p

        path = os.path.join(tmp_path, "no_cutoffs.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)

        cm = ConfigManager(cfg_path=str(path))
        assert cm.config_args['cutoffs'] == {}

    def test_check_missing_keys_raises(self, tmp_path, get_clean_cfg):
        import configparser
        cfg_dict = get_clean_cfg
        del cfg_dict['generals']['n_procs']

        cfg_obj = configparser.ConfigParser()
        for s, p in cfg_dict.items():
            cfg_obj.add_section(s)
            for k, v in p.items():
                cfg_obj[s][k] = str(v)

        path = os.path.join(tmp_path, "missing_key.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)

        with pytest.raises(KeyError, match="is missing"):
            ConfigManager(cfg_path=path)

    def test_read_annotations_paths(self, tmp_path, get_clean_cfg,
                                    write_modified_cfg):
        # --- 1. Absolute Path Test ---
        cfg_abs = get_clean_cfg
        cfg_abs['topo-traj']['topology'] = __file__
        cfg_abs['topo-traj']['trajectory'] = "dummy"
        cfg_abs['generals']['output_dir'] = os.path.join(tmp_path,
                                                         "out_ann_abs")  # UNIQUE DIR 1

        af = os.path.join(tmp_path, "a.txt")
        with open(af, "w") as f: f.write("data")
        cfg_abs['interactions']['annotations'] = str(af)

        path_abs = write_modified_cfg(tmp_path, cfg_abs, 'generals',
                                      'job_name', 'a')
        assert ConfigManager(cfg_path=path_abs).config_args[
                   'annotations'] == str(af)

        # --- 2. Relative Path Test ---
        cfg_rel = get_clean_cfg  # Grab a fresh copy just to be safe
        cfg_rel['topo-traj']['topology'] = __file__
        cfg_rel['topo-traj']['trajectory'] = "dummy"
        cfg_rel['generals']['output_dir'] = os.path.join(tmp_path,
                                                         "out_ann_rel")  # UNIQUE DIR 2

        rf = "r.txt"
        with open(os.path.join(tmp_path, rf), "w") as f: f.write("data")
        cfg_rel['interactions']['annotations'] = rf

        path_rel = write_modified_cfg(tmp_path, cfg_rel, 'generals',
                                      'job_name', 'r')
        assert ConfigManager(cfg_path=path_rel).config_args[
            'annotations'].endswith(rf)

    def test_read_annotations_final_error(self, tmp_path, get_clean_cfg,
                                          write_modified_cfg):
        cfg_dict = get_clean_cfg
        cfg_dict['topo-traj']['topology'] = __file__
        cfg_dict['topo-traj']['trajectory'] = "dummy"
        cfg_dict['interactions']['annotations'] = "completely_fake_path.txt"
        cfg_dict['generals']['output_dir'] = os.path.join(tmp_path,
                                                          "out_ann_err")

        path = write_modified_cfg(tmp_path, cfg_dict, 'generals', 'job_name',
                                  'ann_err')
        with pytest.raises(ValueError, match="Error parsing the key"):
            ConfigManager(cfg_path=path)
