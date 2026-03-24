"""
Creator: rglez
Date: 3/22/26
Env:
Desc:
"""
import configparser
import os
from os.path import dirname, exists, join, normpath

import numba
import pytest

from intermap.managers.config import ConfigManager


class TestDirectoryAndPaths:
    """Tests involving physical paths, directories, and side effects on disk."""

    def test_out_dir_hierarchy(self, tmp_path, get_clean_cfg,
                               write_modified_cfg):
        cfg = get_clean_cfg
        path = write_modified_cfg(tmp_path, cfg, 'generals', 'output_dir',
                                  "results")

        cm = ConfigManager(cfg_path=path)
        assert cm.config_args['output_dir'] == normpath(
            join(dirname(path), "results"))

    def test_dir_creation(self, tmp_path, get_clean_cfg, write_modified_cfg):
        cfg = get_clean_cfg
        out_dir = join(tmp_path, "new_out")
        cfg['generals']['output_dir'] = out_dir

        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name',
                                  'test_job')
        ConfigManager(cfg_path=path)
        assert exists(out_dir)
        assert exists(join(out_dir, "test_job_InterMap.cfg"))

    def test_dir_already_exists(self, tmp_path, get_clean_cfg,
                                write_modified_cfg):
        cfg = get_clean_cfg
        out_dir = join(tmp_path, "already_here")
        os.makedirs(out_dir, exist_ok=True)
        cfg['generals']['output_dir'] = out_dir

        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name',
                                  'test')
        with pytest.raises(FileExistsError):
            ConfigManager(cfg_path=path)

    def test_read_annotations_paths(self, tmp_path, get_clean_cfg,
                                    write_modified_cfg):
        # 1. Absolute Path
        cfg_abs = get_clean_cfg
        cfg_abs['generals']['output_dir'] = join(tmp_path, "out_ann_abs")

        af = join(tmp_path, "a.txt")
        with open(af, "w") as f:
            f.write("data")
        cfg_abs['interactions']['annotations'] = str(af)

        path_abs = write_modified_cfg(tmp_path, cfg_abs, 'generals',
                                      'job_name', 'a')
        assert ConfigManager(cfg_path=path_abs).config_args[
                   'annotations'] == str(af)

        # 2. Relative Path
        cfg_rel = get_clean_cfg
        cfg_rel['generals']['output_dir'] = join(tmp_path, "out_ann_rel")

        rf = "r.txt"
        with open(join(tmp_path, rf), "w") as f:
            f.write("data")
        cfg_rel['interactions']['annotations'] = rf

        path_rel = write_modified_cfg(tmp_path, cfg_rel, 'generals',
                                      'job_name', 'r')
        assert ConfigManager(cfg_path=path_rel).config_args[
            'annotations'].endswith(rf)


class TestStructuralConstraints:
    """
    Tests involving missing sections, missing keys, or structural file errors.
    """

    def test_missing_section(self, tmp_path, get_clean_cfg):
        cfg_obj = configparser.ConfigParser()
        for s, params in get_clean_cfg.items():
            if s != 'generals':
                cfg_obj[s] = params

        path = join(tmp_path, "missing_sec.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)

        with pytest.raises(ValueError, match="Missing sections"):
            ConfigManager(cfg_path=path)

    def test_missing_cutoffs_section(self, tmp_path, get_clean_cfg):
        cfg_dict = get_clean_cfg
        cfg_obj = configparser.ConfigParser()
        cfg_obj.optionxform = str
        for s, p in cfg_dict.items():
            if s != 'cutoffs':
                cfg_obj[s] = p

        path = join(tmp_path, "no_cutoffs.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)

        cm = ConfigManager(cfg_path=str(path))
        assert cm.config_args['cutoffs'] == {}

    def test_check_missing_keys_raises(self, tmp_path, get_clean_cfg):
        cfg_dict = get_clean_cfg
        del cfg_dict['generals']['n_procs']

        cfg_obj = configparser.ConfigParser()
        for s, p in cfg_dict.items():
            cfg_obj.add_section(s)
            for k, v in p.items():
                cfg_obj[s][k] = str(v)

        path = join(tmp_path, "missing_key.cfg")
        with open(path, "w") as f:
            cfg_obj.write(f)

        with pytest.raises(KeyError, match="is missing"):
            ConfigManager(cfg_path=path)

    def test_read_annotations_final_error(self, tmp_path, get_clean_cfg,
                                          write_modified_cfg):
        cfg_dict = get_clean_cfg
        cfg_dict['interactions']['annotations'] = "completely_fake_path.txt"
        cfg_dict['generals']['output_dir'] = join(tmp_path, "out_ann_err")

        path = write_modified_cfg(tmp_path, cfg_dict, 'generals', 'job_name',
                                  'ann_err')
        with pytest.raises(ValueError, match="Error parsing the key"):
            ConfigManager(cfg_path=path)

    def test_read_annotations_true_error(self, tmp_path, get_clean_cfg,
                                         write_modified_cfg):
        cfg_dict = get_clean_cfg
        cfg_dict['interactions']['annotations'] = "True"
        cfg_dict['generals']['output_dir'] = join(tmp_path, "out_ann_true_err")

        path = write_modified_cfg(tmp_path, cfg_dict, 'generals', 'job_name',
                                  'ann_true_err')
        with pytest.raises(ValueError,
                           match="must be set either to False or to a valid path"):
            ConfigManager(cfg_path=path)


class TestMasterParameterValidation:
    """Parametrized tests for all static value boundaries, types, and choices."""

    @pytest.mark.parametrize("section, key, value, expected_error", [
        # --- generals
        ("generals", "n_procs", "0", "out of"),
        ("generals", "n_procs", "abc", "invalid literal for int"),
        ("generals", "n_samples", "0", "out of"),
        ("generals", "n_samples", "abc", "invalid literal for int"),
        ("generals", "n_factor", "-0.1", "out of"),
        ("generals", "n_factor", "abc", "could not convert string to float"),

        # --- topo-traj
        ("topo-traj", "start", "-1", "out of"),
        ("topo-traj", "start", "1.5", "invalid literal for int"),
        ("topo-traj", "stride", "0", "out of"),
        ("topo-traj", "stride", "abc", "invalid literal for int"),
        ("topo-traj", "last", "-2", "out of"),
        ("topo-traj", "last", "any", "invalid literal for int"),
        ("topo-traj", "chunk_size", "0", "out of"),
        ("topo-traj", "chunk_size", "0.5", "invalid literal for int"),

        # --- interactions
        ("interactions", "min_prevalence", "100.1", "out of"),
        ("interactions", "min_prevalence", "-0.1", "out of"),
        ("interactions", "min_prevalence", "none",
         "could not convert string to float"),
        ("interactions", "resolution", "residue_set", "resolution"),
        ("interactions", "interactions", "MagicBond", "Invalid interaction"),

        # --- cutoffs
        ("cutoffs", "fake_cutoff", "5.0", "is not a valid cutoff name"),
    ])
    def test_parameter_boundaries(self, tmp_path, get_clean_cfg,
                                  write_modified_cfg,
                                  section, key, value, expected_error):
        cfg = get_clean_cfg
        # Unique output dir for every case to avoid FileExistsError
        cfg['generals']['output_dir'] = join(tmp_path, f"out_{section}_{key}")
        path = write_modified_cfg(tmp_path, cfg, section, key, value)

        with pytest.raises(ValueError, match=expected_error):
            ConfigManager(cfg_path=path)

    def test_path_parameter_validation(self, tmp_path, get_clean_cfg,
                                       write_modified_cfg):
        """
        Specifically tests PathParam checking physical existence.
        """
        cfg = get_clean_cfg
        cfg['topo-traj']['topology'] = "/non/existent/structure.pdb"

        path = write_modified_cfg(tmp_path, cfg, 'topo-traj', 'topology',
                                  cfg['topo-traj']['topology'])

        with pytest.raises(ValueError, match="No such file"):
            ConfigManager(cfg_path=path)


class TestRuntimeEnforcement:
    """
    Verifies that config parameters actually change the global library states.
    """

    def test_numba_threads_respected(self, tmp_path, get_clean_cfg,
                                     write_modified_cfg):
        """
        Regression test: Ensures n_procs is passed to numba.set_num_threads.
        """
        # 1. Save original state
        original_threads = numba.get_num_threads()

        try:
            # 2. Setup config with a restricted number of threads (1)
            cfg = get_clean_cfg
            target_threads = 5
            path = write_modified_cfg(tmp_path, cfg, 'generals', 'n_procs',
                                      str(target_threads))

            # 3. Initialize Manager (this should trigger numba.set_num_threads)
            ConfigManager(cfg_path=path)

            # 4. Assert Numba global state was updated
            assert numba.get_num_threads() == target_threads

        finally:
            # 5. Always restore the original count for the rest of the tests
            numba.set_num_threads(original_threads)

    def test_dir_creation(self, tmp_path, get_clean_cfg, write_modified_cfg):
        cfg = get_clean_cfg
        out_dir = join(tmp_path, "new_out")
        cfg['generals']['output_dir'] = out_dir

        # Define a specific job name to test naming logic
        my_custom_job = "exp_v1_run_01"
        path = write_modified_cfg(tmp_path, cfg, 'generals', 'job_name',
                                  my_custom_job)

        ConfigManager(cfg_path=path)

        # Assert the filename correctly uses the job_name as a prefix
        expected_filename = f"{my_custom_job}_InterMap.cfg"
        assert exists(out_dir)
        assert exists(join(out_dir, expected_filename))
