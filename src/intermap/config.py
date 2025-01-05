# Created by roy.gonzalez-aleman at 13/11/2023
import configparser
import os
import sys
from collections import Counter
from os.path import abspath, dirname, isabs, join, normpath

import intermap.commons as cmn

inf_int = sys.maxsize
inf_float = float(inf_int)

#: Allowed section templates in the config file
allowed_templates = {'standard':
                         {'generals', 'topo-traj', 'interactions', 'cutoffs'}}

#:  Allowed keys in the config file (dtypes & expected values)
allowed_parameters = {
    # ____ generals
    'generals': {
        'output_dir': {'dtype': 'path', 'check_exist': False},
        'n_procs': {'dtype': int, 'min': 1, 'max': inf_int},
        'job_name': {'dtype': 'path', 'check_exist': False}},
    # ____ topo-traj
    'topo-traj': {
        'topology': {'dtype': 'path', 'check_exist': True},
        'trajectory': {'dtype': 'path', 'check_exist': True},
        'start': {'dtype': int, 'min': 0, 'max': inf_int},
        'last': {'dtype': int, 'min': 0, 'max': inf_int},
        'stride': {'dtype': int, 'min': 1, 'max': inf_int},
        'chunk_size': {'dtype': int, 'min': 1, 'max': inf_int}},
    # ____ interactions
    'interactions': {
        'selection_1': {'dtype': str},
        'selection_2': {'dtype': str},
        'interactions': {'dtype': str,
                         'values': {'all', 'CloseContacts', 'VdWContact',
                                    'Hydrophobic', 'Anionic', 'Cationic',
                                    'MetalDonor', 'MetalAcceptor',
                                    'HBAcceptor',
                                    'HBDonor', 'XBAcceptor', 'XBDonor',
                                    'PiStacking', 'FaceToFace', 'EdgeToFace',
                                    'PiCation', 'CationPi'}}},
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

    def __init__(self, config_path, legal_params, legal_templates):

        # Parse class args
        self.config_path = cmn.check_path(config_path)
        self.legal_params = legal_params
        self.legal_templates = legal_templates

        # Parsing from class args
        self.config_dir = abspath(dirname(self.config_path))
        self.keyless_sections = self.detect_keyless_sections()
        self.config_obj = self.read_config_file()
        self.template = self.detect_template()

        # Run checks
        self.check_missing_keys()
        self.config_args = self.check_params()
        self.parse_and_check_constraints()
        self.config_args['template'] = self.template
        self.check_hydrogens_presence()

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

    def detect_template(self):
        """
        Detect the template of the configuration file

        Raises:
            ValueError: if the workflow or sections are not valid

        Returns:
            workflow: the workflow specified in the configuration
        """
        # Check if the workflow is valid
        try:
            workflow = self.config_obj['generals']['workflow']
        except KeyError:
            raise KeyError(
                'The configuration file must have a [generals] section with'
                ' a "workflow" key specifying the workflow to run.')

        legal_workflows = self.legal_templates.keys()
        if workflow not in legal_workflows:
            raise ValueError(
                f'Workflow "{workflow}" is not a valid one. '
                f'Currently allowed workflows are: {legal_workflows}')

        # Check if the sections are valid
        valid_sections = self.legal_templates[workflow]
        current_sections = set(self.config_obj.sections())

        # Check if the sections are valid for the workflow specified
        diff1 = set.difference(valid_sections, current_sections)
        if diff1:
            raise ValueError(
                f'The configuration file is not valid for the workflow'
                f' "{workflow}". The following sections are missing: {diff1}')

        diff2 = set.difference(current_sections, valid_sections)
        if diff2:
            raise ValueError(
                f'The configuration file is not valid for the workflow'
                f' "{workflow}". The following sections are unnecesary: {diff2}')

        return workflow

    def check_missing_keys(self):
        """
        Check for missing keys in the configuration file

        Raises:
            KeyError: if a key is missing in the configuration file
        """
        current_template = self.legal_templates[self.template].copy()
        current_params = self.legal_params
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

                param_info = self.legal_params[section][key]
                dtype = param_info['dtype']
                if dtype in {float, int}:
                    param_obj = NumericParam(key, dtype(value), **param_info)
                elif dtype == 'path':
                    value = value if isabs(value) else join(config_dir, value)
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


class STDConfig(Config):
    """
    Specific parser for STDock's config files. It inherits from a more general
    config parser and then perform STDock-related checkings.
    """

    def parse_and_check_constraints(self):
        config_sections = self.config_obj.sections()

        # 1. todo: Check [docking] section

        # 2. Check [std-spectra] sectioan
        if 'std-spectra' in config_sections:
            self.config_args['std-spectra'] = self.parse_spectra()

        # 3. Check [std-regions] section
        if 'std-regions' in config_sections:
            self.config_args['std-regions'] = self.parse_regions()

        # 4. Check [std-epitopes] section
        if 'std-epitope' in config_sections:
            self.config_args['std-epitopes'] = self.parse_epitopes()

        # 1. Build dir hierarchy
        self.build_dir_hierarchy()

    def build_dir_hierarchy(self):
        """
        Build STDock directory hierarchy
        """
        # If output_dir exists, raise
        outdir = self.config_args['output_dir']

        try:
            os.makedirs(outdir)
        except FileExistsError:
            raise FileExistsError(
                f'The output directory {outdir} already exists. Please, '
                f'choose another one, or delete the existing one.')

        for dir_name in ['STD', 'DOCKING']:
            self.config_args[dir_name] = join(outdir, dir_name)
            os.makedirs(self.config_args[dir_name])

        # Write the configuration file for reproducibility
        stdock_config = join(self.config_args['output_dir'], 'stdock-job.cfg')
        with open(stdock_config, 'wt') as ini:
            self.config_obj.write(ini)

    def parse_spectra(self):
        """
        Parse the [std-spectra] section

        Returns:
            spectra_dict: a nested dict with ligand concentrations and
                          saturation times respectively
        """
        # Extract keys and values
        config = self.config_obj
        keys = list(config['std-spectra'].keys())
        spectra = [x.strip() for x in list(config['std-spectra'].values())]

        # Perform path validity checks for spectra files
        spectra_paths = [cmn.check_path(x) for x in spectra]

        # Get types, times & concentrations in the [std-spectra] section
        types, times, concs = [], [], []
        for key in keys:
            splitted = key.split('_')
            if len(splitted) != 3:
                raise ValueError(
                    'Keys of the [std-spectra] section of the config file must'
                    ' be labeled as [on, off or diff]_[ligand_concentration]_[saturation_time].'
                    ' For example: on_0.1_0.5, off_0.1_0.5, diff_0.1_0.5')

            types.append(splitted[0])
            times.append(float(splitted[1]))
            concs.append(float(splitted[2]))

        # Check types of spectra are: on, off or diff
        types_counts = Counter(types)
        valid_types = ['on', 'off', 'diff']
        types_are_valid = all([x in valid_types for x in types_counts.keys()])
        if not types_are_valid:
            raise ValueError(
                f'Invalid spectrum types in the [std-spectra] of the '
                f'config file. Valid types are: {valid_types}')

        # Check only two types are present
        if len(types_counts) != 2:
            raise ValueError(
                'There must be only two types of spectrum specified in the '
                '[std-spectra] section of the config file: (on and off) OR'
                ' (diff and off)')

        # Check types number consistency
        if len(set(types_counts.values())) != 1:
            raise ValueError(
                'There must be the same number of (on and off) OR (diff and off)'
                ' spectra specified in the [std-spectra] section of the config'
                ' file')

        # Check types consistency
        conc_times = Counter(list(zip(concs, times))).values()
        if set(conc_times) != {2}:
            raise ValueError(
                'There must be two spectrum per each'
                ' (concentration/saturation_time) level')

        # Check times consistency
        if len(set(times)) < 2:
            raise ValueError(
                'There must be at least two different saturation times for each '
                'ligand concentration specified')

        # Check concs consistency
        if len(set(concs)) < 1:
            raise ValueError(
                'There must be at least one ligand concentration specified')

        if len(set(Counter(concs).values())) != 1:
            raise ValueError(
                'There must be the same number of spectra for each ligand '
                'concentration specified in the [std-spectra] section of the '
                'config file')

        # Get all the spectra organized as a dict
        spectra_dict = cmn.recursive_defaultdict()
        for i, key in enumerate(keys):
            splitted = key.split('_')
            type_ = splitted[0]
            lig_conc, sat_time = map(float, splitted[2:0:-1])
            spectra_dict[lig_conc][sat_time][type_] = spectra_paths[i]
        return spectra_dict

    def parse_regions(self):
        """
        Parse the [std-regions] section

        Returns:
            regions_dict: a dict with regions limits to integrate
        """

        # Extract keys and values
        config = self.config_obj
        keys_raw = list(config['std-regions'].keys())
        keys = [cmn.check_epitope_label(x) for x in keys_raw]
        values_raw = list(config['std-regions'].values())

        # Get all the integral regions organized as a dict
        regions_dict = {}
        for i, key in enumerate(keys):
            upper_raw, lower_raw = values_raw[i].split(',')
            if (upper := float(upper_raw)) > (lower := float(lower_raw)):
                regions_dict.update(
                    {i: {'label': key, 'upper': upper, 'lower': lower}})
            else:
                raise ValueError(
                    f'Upper limit must be greater than the lower limit in the'
                    f' {key} region.')
        return regions_dict

    def parse_epitopes(self):
        """
        Parse the [std-epitope] section

        Returns:
            epitope: a dict with the epitope values
        """
        config = self.config_obj
        keys_raw = list(config['std-epitope'].keys())
        keys = [cmn.check_epitope_label(x) for x in keys_raw]
        values_raw = [float(x) for x in config['std-epitope'].values()]
        epitope = dict(zip(keys, values_raw))
        return epitope

# %%===========================================================================
# Debugging area
# =============================================================================
# configs = {
#     # 'map-from-spectra':
#     #     '/home/gonzalezroy/RoyHub/stdock/tests/paper/trolls/config.cfg',
#     # 'map-from-spectra-then-dock-small':
#     #     '/home/gonzalezroy/RoyHub/stdock/tests/paper/imaging-docking/config_kd_docking.cfg',
#     'map-from-spectra-then-dock-small':
#         '/home/gonzalezroy/RoyHub/stdock/tests/paper/imaging-docking/config_kd_docking.cfg',
# }
#
#
# def check_workflow(workflow, configs):
#     """
#     Check a specific workflow
#
#     Args:
#         workflow: stdock workflow
#         configs: dict with the config files
#
#     Returns:
#
#     """
#     params = allowed_parameters
#     templates = allowed_templates
#     config_path = configs[workflow]
#     self = STDConfig(config_path, params, templates)
#     args = self.config_args
#     return args
#
#
# workflow = 'map-from-spectra-then-dock-small'
# args = check_workflow(workflow, configs)
