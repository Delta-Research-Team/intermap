# Created by roy.gonzalez-aleman at 13/11/2023
import configparser
import logging
import os
import re
import sys
from os.path import abspath, basename, dirname, isabs, join, normpath

import numpy as np

import intermap.commons as cmn
import intermap.managers.cutoffs as cf

inf_int = sys.maxsize
inf_float = float(inf_int)

proj_dir = os.sep.join(dirname(os.path.abspath(__file__)).split(os.sep)[:-2])


# =============================================================================
# Allowed sections & parameters
# =============================================================================
def print_colored_ascii():
    """
    Print the InterMap logo in colored ASCII art.
    """
    html_path = join(proj_dir, 'intermap', 'binary_imap.html')

    try:
        with open(html_path, encoding='utf-8') as file:
            content = file.read()
            pattern = r'<div style="margin: 20px 0;">(.*?)</div>'
            match = re.search(pattern, content, re.DOTALL)

            if match:
                art_content = match.group(1)
                art_content = art_content.replace('<br/>', '\n').replace(
                    '<br>', '\n')
                art_content = re.sub(
                    r'<span style="color: rgb\((\d+),\s*(\d+),\s*(\d+)\)">(.*?)</span>',
                    lambda
                        m: f"\033[38;2;{m.group(1)};{m.group(2)};{m.group(3)}m{m.group(4)}\033[0m",
                    art_content
                )

                art_content = re.sub(r'<[^>]+>', '', art_content)
                art_content = art_content.replace('&nbsp;', ' ')
                lines = art_content.split('\n')
                formatted_lines = []
                for line in lines:
                    if line.strip():
                        if not line.endswith('\033[0m'):
                            line += '\033[0m'
                        formatted_lines.append(line)

                print('\n\n')
                print('\n'.join(formatted_lines))
                print('\033[0m')

    except Exception as e:
        print(f"Error when reading HTML: {e}")


def detect_config_path(mode='debug'):
    """
    Detect the configuration file path

    Args:
        mode: running mode

    Returns:
        config_path: path to the configuration file
    """
    if mode == 'production':
        if len(sys.argv) != 2:
            raise ValueError(
                '\nInterMap syntax is: intermap path-to-config-file')
        config_path = sys.argv[1]
    elif mode == 'debug':
        # config_path = '/media/gonzalezroy/Expansion/romie/TRAJECTORIES_INPUTS_DATA_mpro_wt_variants_amarolab/a173v/imap.cfg'
        config_path = '/media/rglez/Roy2TB/Dropbox/DockingBoys/storage/identity_prolif_100frames/IgG3/imap.cfg'
    else:
        raise ValueError('Only modes allowed are production and running')
    return config_path


def start_logger(log_path):
    """
    Start the logger for the InterMap run.

    Args:
        log_path (str): Path to the log file.

    Returns:
        logger (logging.Logger): Logger object.
    """
    logger = logging.getLogger('InterMapLogger')
    logger.setLevel("DEBUG")

    console_handler = logging.StreamHandler()
    console_handler.setLevel("DEBUG")
    formatter = logging.Formatter(
        ">>>>>>>> {asctime} - {levelname} - {message}\n",
        style="{",
        datefmt="%Y-%m-%d %H:%M")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_path, encoding="utf-8", mode='w')
    file_handler.setLevel("DEBUG")
    logger.addHandler(file_handler)
    return logger


class Param:
    """
    Base class for parameter checking
    """

    def __init__(self, key, value, *args, **kwargs):
        self.key = key
        self.value = value
        self.args = args
        self.kwargs = kwargs

    def check(self):
        """
        Check the parameter
        """
        raise NotImplementedError()


class NumericParam(Param):
    """
    Check numeric parameters
    """

    def check(self):
        dtype = self.kwargs['dtype']
        minim = self.kwargs['min']
        maxim = self.kwargs['max']
        cmn.check_numeric_in_range(self.key, self.value, dtype, minim, maxim)


class PathParam(Param):
    """
    Check path
    """

    def check(self):
        path = self.value
        cmn.check_path(path, check_exist=self.kwargs['check_exist'])


class ChoiceParam(Param):
    """
    Check choices
    """

    def check(self):
        choices = self.kwargs['values']
        if choices and (not (self.value in choices)):
            raise ValueError(
                f'\n Error in {self.key}. Passed "{self.value}" but available'
                f' options are: {choices}.')


class Config:
    """
    Base class for config file parsing
    """
    allowed_parameters = {
        # ____ generals
        'generals': {
            'output_dir': {'dtype': 'path', 'check_exist': False},
            'n_procs': {'dtype': int, 'min': 1, 'max': inf_int},
            'job_name': {'dtype': 'path', 'check_exist': False},
            'n_samples': {'dtype': int, 'min': 1, 'max': inf_int},
            'n_factor': {'dtype': float, 'min': 1, 'max': inf_float},
        },

        # ____ topo-traj
        'topo-traj': {
            'topology': {'dtype': 'path', 'check_exist': True},
            'trajectory': {'dtype': str, 'values': None},
            'start': {'dtype': int, 'min': 0, 'max': inf_int},
            'last': {'dtype': int, 'min': -1, 'max': inf_int},
            'stride': {'dtype': int, 'min': 1, 'max': inf_int},
            'chunk_size': {'dtype': int, 'min': 1, 'max': inf_int}},

        # ____ interactions
        'interactions': {
            'selection_1': {'dtype': str, 'values': None},
            'selection_2': {'dtype': str, 'values': None},
            'min_prevalence': {'dtype': float, 'min': 0, 'max': 100},
            'interactions': {'dtype': str, 'values': None},
            'export_csv': {'dtype': str, 'values': {'True', 'False'}},
            'resolution': {'dtype': str, 'values': {'atom', 'residue'}},
            'annotations': {'dtype': str, 'values': None},
            'format': {'dtype': str, 'values': {'simple', 'extended'}}},

        # ____ cutoffs
        'cutoffs': None
    }

    def __init__(self, mode='production'):

        # Print the header
        print_colored_ascii()

        # Detect config
        self.config_path = cmn.check_path(detect_config_path(mode=mode))
        self.legal_params = self.allowed_parameters

        # Parsing from class args
        self.config_dir = abspath(dirname(self.config_path))
        self.keyless_sections = self.detect_keyless_sections()
        self.config_obj = self.read_config_file()

        # Run checks
        self.check_missing_keys()
        self.config_args = self.check_params()
        self.parse_and_check_constraints()

    def detect_keyless_sections(self):
        """
        Detect sections without keys in the configuration file

        Returns:
            keyless_sections: a list with the sections without keys
        """
        params = self.legal_params
        keyless_sections = [x for x in params if params[x] is None]
        return keyless_sections

    def read_config_file(self):
        """
        Read the configuration file

        Returns:
            config_obj: a ConfigParser object of the configuration file
        """
        config_obj = configparser.ConfigParser(allow_no_value=True,
                                               inline_comment_prefixes='#')
        config_obj.optionxform = str
        config_obj.read(self.config_path)
        return config_obj

    def check_missing_keys(self):
        """
        Check for missing keys in the configuration file

        Raises:
            KeyError: if a key is missing in the configuration file
        """
        current_params = self.legal_params
        current_template = list(current_params.keys())

        [current_template.remove(x) for x in self.keyless_sections if
         x in current_template]

        for section in current_template:
            config_file_keys = list(self.config_obj[section].keys())
            for key in current_params[section]:
                if key not in config_file_keys:
                    raise KeyError(
                        f'Key "{key}" is missing from the "{section}" section'
                        f' of the configuration file. Please set its value.')

    def check_params(self):
        """
        Check the parameters in the configuration file

        Returns:
            config_args: a dict with the parsed and checked parameters
        """
        config_args = dict()

        parsed_sections = self.config_obj.sections().copy()
        [parsed_sections.remove(x) for x in self.keyless_sections if
         x in parsed_sections]

        for section in parsed_sections:
            items = self.config_obj[section].items()
            for key, value in items:

                try:
                    param_info = self.legal_params[section][key]
                except KeyError:
                    raise KeyError(
                        f'Key "{key}" is forbidden in the section "{section}".')

                dtype = param_info['dtype']
                if dtype in {float, int}:
                    param_obj = NumericParam(key, dtype(value), **param_info)
                elif dtype == 'path':
                    absolute = normpath(join(dirname(self.config_path), value))
                    value = value if isabs(value) else absolute
                    param_obj = PathParam(key, value, **param_info)
                elif dtype == str:
                    param_obj = ChoiceParam(key, value, **param_info)
                else:
                    raise ValueError(
                        f"\n{section}.{key}'s dtype is wrong: {dtype}.")
                param_obj.check()

                parsed_value = normpath(value) if dtype == 'path' else dtype(
                    value)
                config_args.update({key: parsed_value})

        return config_args

    def parse_and_check_constraints(self):
        """Check for specific constraints in the STDock config file
        """
        raise NotImplementedError


class ConfigManager(Config):
    """
    Specific parser for STDock's config files. It inherits from a more general
    config parser and then perform STDock-related checkings.
    """

    def parse_and_check_constraints(self):
        # 1. Build the directory hierarchy
        self.build_dir_hierarchy()

        # 2. Parse the cutoffs
        self.parse_cutoffs()

        # 3. Parse the interactions
        self.parse_interactions()

        # 4. Parse the annotations
        self.read_annotations()

        # 5. Start logging
        args = self.config_args
        base_name = basename(args['job_name'])
        log_path = join(args['output_dir'], f"{base_name}_InterMap.log")
        logger = start_logger(log_path)

        value = args['trajectory']
        ntrajs = 'Trajectory' if len(value.split(',')) == 1 else 'Trajectories'

        logger.info(
            f"Starting InterMap with the following static parameters:\n"
            f"\n Job name: {args['job_name']}"
            f"\n Topology: {args['topology']}"
            f"\n {ntrajs}: {value}"
            f"\n String for selection #1: {args['selection_1']}"
            f"\n String for selection #2: {args['selection_2']}"
            f"\n Output directory: {args['output_dir']}"
            f"\n Chunk size: {args['chunk_size']}"
            f"\n Number of processors: {args['n_procs']}"
            f"\n Min prevalence: {args['min_prevalence']}"
            f"\n Report's format: {args['format']}"
            f"\n Resolution: {args['resolution']}\n"
        )

    def build_dir_hierarchy(self):
        """
        Build STDock directory hierarchy
        """
        # If output_dir exists, raise
        outdir = self.config_args['output_dir']

        try:
            os.makedirs(outdir, exist_ok=True)
        except FileExistsError:
            raise FileExistsError(
                f'The output directory {outdir} already exists. Please, '
                f'choose another one, or delete the existing.')

        # Write the configuration file for reproducibility
        config = join(self.config_args['output_dir'], 'InterMap-job.cfg')
        with open(config, 'wt') as ini:
            self.config_obj.write(ini)

    def parse_cutoffs(self):
        """
        Parse the cutoffs
        """

        # Get the internal cutoffs
        internal_names = list(cf.cutoffs.keys())

        # Parse the cutoffs from the config file
        try:
            cutoffs = self.config_obj['cutoffs']
        except KeyError:
            cutoffs = {}

        # Check if the cutoffs are valid
        config_cutoffs = dict()
        for key, value in cutoffs.items():
            if key not in cf.cutoffs:
                raise ValueError(
                    f"{key} is not a valid cutoff name.\n"
                    f"The supported list is:\n"
                    f"{internal_names}")
            config_cutoffs.update({key: float(value)})

        self.config_args.update({'cutoffs': config_cutoffs})

    def parse_interactions(self):
        """
        Parse the interactions from the config file
        """
        raw_inters = self.config_obj['interactions']['interactions']
        if raw_inters == 'all':
            parsed_inters = np.asarray(cf.interactions + ['WaterBridge'])
        else:
            parsed_inters = [x.strip() for x in raw_inters.split(',') if
                             x != '']
            parsed_inters = np.asarray(parsed_inters)

        implemented = cf.interactions + ['WaterBridge']
        for inter in parsed_inters:
            if inter not in implemented:
                raise ValueError(
                    f"Invalid interaction specified: {inter}. The list of"
                    f" currently implemented interactions is:\n {implemented}\n")

        self.config_args.update({'interactions': parsed_inters})

    def read_annotations(self):
        """
        Parse the annotations from the config file
        """
        # Get the absolute path for the annotations file
        raw_value = self.config_obj['interactions']['annotations'].strip()
        cfg_path = self.config_path

        if raw_value == 'False':
            self.config_args.update({'annotations': False})

        else:
            try:
                abs_path = cmn.check_path(raw_value)
                self.config_args.update({'annotations': abs_path})
            except ValueError:
                abs_path = normpath(join(dirname(cfg_path), raw_value))

                try:
                    abs_path = cmn.check_path(abs_path)
                    self.config_args.update({'annotations': abs_path})
                except ValueError:
                    raise ValueError(
                        f'Error parsing the key "annotations" in the '
                        f'"interactions" section of the config file. The value'
                        f' must be set either to False or to a valid path. The'
                        f' provided value was: {raw_value}.\n')

# %%===========================================================================
# Debugging area
# =============================================================================
# config_path = '/home/gonzalezroy/RoyHub/NUC-STRESS-RGA/data/0A-prelude/DrHU-repliques/second_step/imap.cfg'
# self = ConfigManager(mode='debug')
