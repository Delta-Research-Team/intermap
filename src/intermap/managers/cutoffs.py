# Created by rglez at 12/9/24
"""
Cutoff distances for the different definitions of interactions
"""
import logging
import re
from pprint import pformat

import numpy as np
import numpy_indexed as npi
from numba.typed import List

logger = logging.getLogger('InterMapLogger')

cutoffs = {

    'dist_cut_CloseContacts': 3.0,
    'dist_cut_Hydrophobic': 4.5,
    'dist_cut_Ionic': 4.5,
    'dist_cut_Metalic': 2.8,

    'dist_cut_HA': 3.5,
    'dist_cut_DA': 3.5,
    'min_ang_DHA': 130,
    'max_ang_DHA': 180,

    'dist_cut_XA': 3.5,
    'dist_cut_XD': 3.5,
    'min_ang_DXA': 0,
    'max_ang_DXA': 90,

    'dist_cut_FaceToFace': 5.5,
    'min_ang_nn_FaceToFace': 0,
    'max_ang_nn_FaceToFace': 35,
    'min_ang_nc_FaceToFace': 0,
    'max_ang_nc_FaceToFace': 30,

    'dist_cut_EdgeToFace': 6.5,
    'min_ang_nn_EdgeToFace': 50,
    'max_ang_nn_EdgeToFace': 90,
    'min_ang_nc_EdgeToFace': 0,
    'max_ang_nc_EdgeToFace': 30,
    'intersect_radius_EdgeToFace': 1.5,

    'dist_cut_PiStacking': 6.5,
    'min_ang_nn_PiStacking': 0,
    'max_ang_nn_PiStacking': 90,
    'min_ang_nc_PiStacking': 0,
    'max_ang_nc_PiStacking': 30,

    'dist_cut_PiCation': 4.5,
    'min_ang_PiCation': 0,
    'max_ang_PiCation': 30,
}

interactions = [
    'CloseContacts', 'VdWContact', 'Hydrophobic', 'Anionic',
    'Cationic', 'MetalDonor', 'MetalAcceptor', 'HBAcceptor',
    'HBDonor', 'XBAcceptor', 'XBDonor', 'FaceToFace', 'EdgeToFace',
    'PiStacking', 'PiCation', 'CationPi']


def get_cutoff(cutoff_name, args=()):
    """
    Parse the cutoff name from args or from internal values if not found in args.

    Args:
        args (Namespace): Arguments from the parser.
        cutoff_name (str): Name of the cutoff.

    Returns:
        float: Value of the cutoff.
    """

    # Raise an error if the cutoff name is not in the class
    if cutoff_name not in cutoffs:
        raise ValueError(
            f"{cutoff_name} is not a valid cutoff name.\n"
            f" The supported list is:\n"
            f"{[x for x in cutoffs.keys() if not x.startswith('__')]}")

    # Get the value from the config file
    elif cutoff_name in cutoffs:
        return args[cutoff_name]

    # Get the value from the class
    else:
        return cutoffs[cutoff_name]


class CutoffsManager:
    """
    Class to manage the cutoffs for the different interactions
    """

    def __init__(self, args, iman):
        # Parse args from the config file
        self.args = args
        self.iman = iman

        # Get all the interactions and their parsed cutoffs
        self.max_vdw = iman.get_max_vdw_dist()
        self.parsed_cutoffs = self.parse_cutoffs()
        self.inters_internal, self.cuts_internal = self.get_inters_cutoffs()

        # Split the cutoffs and interactions to compute
        (self.selected_aro, self.selected_others, self.cuts_aro,
         self.cuts_others, self.len_others,
         self.len_aro) = self.split_cutoffs()

        self.max_dist_aro = self.cuts_aro[
                            :2].max() if 'None' not in self.selected_aro else 0

        self.max_dist_others = max(self.cuts_others[:2].max(),
                                   self.max_vdw) if 'None' not in self.selected_others else 0

    def parse_cutoffs(self):
        """
        Parse the cutoffs from the args or the cf module if not found in args.

        Returns:
            dict: Dictionary with the cutoffs.
        """
        cuts = self.args.cutoffs
        cutoffs_raw = {
            # Undirected 1D interactions
            'CloseContacts':
                {'distCut1': get_cutoff('dist_cut_CloseContacts', cuts)},
            'VdWContact':
                {'distCut1': 0},
            'Hydrophobic':
                {'distCut1': get_cutoff('dist_cut_Hydrophobic', cuts)},

            # Directed 1D interactions
            'Anionic':
                {'distCut1': get_cutoff('dist_cut_Ionic', cuts)},
            'Cationic':
                {'distCut1': get_cutoff('dist_cut_Ionic', cuts)},
            'MetalDonor':
                {'distCut1': get_cutoff('dist_cut_Metalic', cuts)},
            'MetalAcceptor':
                {'distCut1': get_cutoff('dist_cut_Metalic', cuts)},

            # Directed 2D1A interactions
            'HBAcceptor':
                {'distCut1': get_cutoff('dist_cut_DA', cuts),
                 'distCut2': get_cutoff('dist_cut_HA', cuts),
                 'minAng1': get_cutoff('min_ang_DHA', cuts),
                 'maxAng1': get_cutoff('max_ang_DHA', cuts)},
            'HBDonor':
                {'distCut1': get_cutoff('dist_cut_DA', cuts),
                 'distCut2': get_cutoff('dist_cut_HA', cuts),
                 'minAng1': get_cutoff('min_ang_DHA', cuts),
                 'maxAng1': get_cutoff('max_ang_DHA', cuts)},
            'XBAcceptor':
                {'distCut1': get_cutoff('dist_cut_XA', cuts),
                 'distCut2': get_cutoff('dist_cut_XD', cuts),
                 'minAng1': get_cutoff('min_ang_DXA', cuts),
                 'maxAng1': get_cutoff('max_ang_DXA', cuts)},
            'XBDonor':
                {'distCut1': get_cutoff('dist_cut_XA', cuts),
                 'distCut2': get_cutoff('dist_cut_XD', cuts),
                 'minAng1': get_cutoff('min_ang_DXA', cuts),
                 'maxAng1': get_cutoff('max_ang_DXA', cuts)},

            # Undirected XD2A interactions
            'FaceToFace':
                {'distCut1': get_cutoff('dist_cut_FaceToFace', cuts),
                 'minAng1': get_cutoff('min_ang_nn_FaceToFace', cuts),
                 'maxAng1': get_cutoff('max_ang_nn_FaceToFace', cuts),
                 'minAng2': get_cutoff('min_ang_nc_FaceToFace', cuts),
                 'maxAng2': get_cutoff('max_ang_nc_FaceToFace', cuts)},
            'EdgeToFace':
                {'distCut1': get_cutoff('dist_cut_EdgeToFace', cuts),
                 'distCut2': get_cutoff('intersect_radius_EdgeToFace', cuts),
                 'minAng1': get_cutoff('min_ang_nn_EdgeToFace', cuts),
                 'maxAng1': get_cutoff('max_ang_nn_EdgeToFace', cuts),
                 'minAng2': get_cutoff('min_ang_nc_EdgeToFace', cuts),
                 'maxAng2': get_cutoff('max_ang_nc_EdgeToFace', cuts), },
            'PiStacking':
                {'distCut1': get_cutoff('dist_cut_PiStacking', cuts),
                 'minAng1': get_cutoff('min_ang_nn_PiStacking', cuts),
                 'maxAng1': get_cutoff('max_ang_nn_PiStacking', cuts),
                 'minAng2': get_cutoff('min_ang_nc_PiStacking', cuts),
                 'maxAng2': get_cutoff('max_ang_nc_PiStacking', cuts)},

            # Directed 1D1A interactions
            'PiCation':
                {'distCut1': get_cutoff('dist_cut_PiCation', cuts),
                 'minAng1': get_cutoff('min_ang_PiCation', cuts),
                 'maxAng1': get_cutoff('max_ang_PiCation', cuts)},
            'CationPi':
                {'distCut1': get_cutoff('dist_cut_PiCation', cuts),
                 'minAng1': get_cutoff('min_ang_PiCation', cuts),
                 'maxAng1': get_cutoff('max_ang_PiCation', cuts)}}

        return cutoffs_raw

    def get_inters_cutoffs(self):
        """
        Get the interaction names and cutoffs from the raw dictionary

        Returns:

        """
        # Parse the cutoffs either from the args or the cf module
        cuts_all = self.parsed_cutoffs
        inters_all = sorted(list(cuts_all.keys()))

        # Create the cutoffs matrix
        N = len(inters_all)
        M = len(set([y for x in cuts_all.values() for y in x]))
        cuts_table = np.zeros((M, N), dtype=np.float32)

        # Fill the cutoffs matrix
        for inter in inters_all:
            if 'distCut1' in cuts_all[inter]:
                cut_value = cuts_all[inter]['distCut1']
                cuts_table[0, inters_all.index(inter)] = cut_value
            if 'distCut2' in cuts_all[inter]:
                cut_value = cuts_all[inter]['distCut2']
                cuts_table[1, inters_all.index(inter)] = cut_value
            if 'minAng1' in cuts_all[inter]:
                cut_value = cuts_all[inter]['minAng1']
                cuts_table[2, inters_all.index(inter)] = cut_value
            if 'maxAng1' in cuts_all[inter]:
                cut_value = cuts_all[inter]['maxAng1']
                cuts_table[3, inters_all.index(inter)] = cut_value
            if 'minAng2' in cuts_all[inter]:
                cut_value = cuts_all[inter]['minAng2']
                cuts_table[4, inters_all.index(inter)] = cut_value
            if 'maxAng2' in cuts_all[inter]:
                cut_value = cuts_all[inter]['maxAng2']
                cuts_table[5, inters_all.index(inter)] = cut_value

        return inters_all, cuts_table

    def split_cutoffs(self):
        """
        Split the cutoffs and interactions to compute in two groups: aromatic and non-aromatic

        Returns:
            to_compute_aro (ndarray): Aromatic interactions to compute
            to_compute_others (ndarray): Non-aromatic interactions to compute
            cutoffs_aro (ndarray): Cutoffs for the aromatic interactions
            cutoffs_others (ndarray): Cutoffs for the non-aromatic interactions
        """

        inters_requested = self.iman.inters_requested

        # Parse aromatics
        bit_aro = [y for x in inters_requested if re.search(r'Pi|Face', x) for
                   y in npi.indices(self.inters_internal, [x])]
        if len(bit_aro) == 0:
            selected_aro = List(['None'])
        else:
            selected_aro = List([self.inters_internal[x] for x in bit_aro])

        # Parse non-aromatics
        bit_others = [y for x in inters_requested if
                      not re.search(r'Pi|Face', x) for y
                      in npi.indices(self.inters_internal, [x])]
        if len(bit_others) == 0:
            selected_others = List(['None'])
        else:
            selected_others = List(
                [self.inters_internal[x] for x in bit_others])

        # Get the cutoffs
        cutoffs_aro = self.cuts_internal[:, bit_aro]
        cutoffs_others = self.cuts_internal[:, bit_others]

        # Get the lengths
        len_others = len(
            selected_others) if 'None' not in selected_others else 0
        len_aro = len(selected_aro) if 'None' not in selected_aro else 0

        cutoffs_str = {x: self.parsed_cutoffs[x] for x in self.parsed_cutoffs
                       if x in inters_requested}
        logger.info(
            f"Interactions to compute & their cutoffs:\n {pformat(cutoffs_str)}")
        return (selected_aro, selected_others, cutoffs_aro, cutoffs_others,
                len_others, len_aro)
