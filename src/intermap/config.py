# Created by roy.gonzalez-aleman at 13/11/2023
import configparser
import os
import sys
from os.path import abspath, dirname, isabs, join, normpath

import numpy as np

import intermap.commons as cmn
import intermap.cutoffs as cf

inf_int = sys.maxsize
inf_float = float(inf_int)
proj_dir = os.sep.join(dirname(os.path.abspath(__file__)).split(os.sep)[:-2])
#: Allowed section templates in the config file

#:  Allowed keys in the config file (dtypes & expected values)
allowed_parameters = {
    # ____ generals
    'generals': {
        'output_dir': {'dtype': 'path', 'check_exist': False},
        'n_procs': {'dtype': int, 'min': 1, 'max': inf_int},
        'job_name': {'dtype': 'path', 'check_exist': False},
    },
    # ____ topo-traj
    'topo-traj': {
        'topology': {'dtype': 'path', 'check_exist': True},
        'trajectory': {'dtype': 'path', 'check_exist': True},
        'start': {'dtype': int, 'min': 0, 'max': inf_int},
        'last': {'dtype': int, 'min': -1, 'max': inf_int},
        'stride': {'dtype': int, 'min': 1, 'max': inf_int},
        'chunk_size': {'dtype': int, 'min': 1, 'max': inf_int}},
    # ____ interactions
    'interactions': {
        'selection_1': {'dtype': str, 'values': None},
        'selection_2': {'dtype': str, 'values': None},
        'min_prevalence': {'dtype': float, 'min': 0, 'max': 100},
        'interactions': {'dtype': str,
                         'values': None},
        'format': {'dtype': str, 'values': {'simple', 'extended'}}},

    # ____ cutoffs
    'cutoffs': None
}


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

    def __init__(self, config_path, legal_params):

        # Parse class args
        self.config_path = cmn.check_path(config_path)
        self.legal_params = legal_params

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
                        f'Key "{key}" is missing from the config'
                        f' file. Please specify its value.')

    def check_params(self):
        """
        Check the parameters in the configuration file

        Returns:
            config_args: a dict with the parsed and checked parameters
        """
        config_args = dict()

        config_dir = self.config_dir
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
                        f'Key "{key}" is not avilable in the section "{section}".')

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


class InterMapConfig(Config):
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
        prefixes = ('min', 'max', 'dist', 'ang')
        internal_names = [x for x in dir(cf) if x.startswith(prefixes)]

        cutoffs = self.config_obj['cutoffs']
        config_cutoffs = dict()
        for key, value in cutoffs.items():
            if key not in cf.__dict__:
                raise ValueError(
                    f"{key} is not a valid cutoff name.\n"
                    f"The supported list is:\n"
                    f"{internal_names}")
            config_cutoffs.update({key: float(value)})

        parsed_cutoffs = cf.parse_cutoffs(config_cutoffs)
        self.config_args.update({'cutoffs': parsed_cutoffs})

    def parse_interactions(self):
        raw_inters = self.config_obj['interactions']['interactions']
        if raw_inters != 'all':
            parsed_inters = [x.strip() for x in raw_inters.split(',')]
            parsed_inters = np.asarray(parsed_inters)
            implemented = set(self.config_args['cutoffs'].keys())
            for inter in parsed_inters:
                if inter not in implemented:
                    raise ValueError(
                        f"Invalid interaction specified: {inter}. The list of"
                        f" currently implemented interactions is:\n {implemented}\n")

            self.config_args.update({'interactions': parsed_inters})

# %%===========================================================================
# Debugging area
# =============================================================================
# config_path = '/home/rglez/RoyHub/intermap/tests/imap_lig-prot.cfg'
# self = InterMapConfig(config_path, allowed_parameters)
