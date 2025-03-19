# Created by rglez at 12/9/24
"""
Cutoff distances for the different definitions of interactions
"""
import sys

import numpy as np

# =============================================================================
# Default cutoffs for the different types of interactions
# =============================================================================

# ==== CloseContacts ==========================================================
dist_cut_CloseContacts = 3.0
dist_cut_Ionic = 4.5
dist_cut_Hydrophobic = 4.5
dist_cut_Metalic = 2.8

# ==== HBonds =================================================================
dist_cut_HA = 2.5
dist_cut_DA = 3.9
min_ang_DHA = 130
max_ang_DHA = 180

# ==== XBonds =================================================================
dist_cut_XA = 2.5
dist_cut_XD = 3.9
ang_cut_DXA = 90

# ==== PiCation ===============================================================
dist_cut_PiCation = 4.5
min_ang_PiCation = 0
max_ang_PiCation = 30

# ==== PiStacking =============================================================
dist_cut_PiStacking = 6.0
min_dist_PiStacking = 3.8
min_ang_PiStacking = 0
max_ang_PiStacking = 90

# ==== EdgeToFace =============================================================
dist_cut_EdgeToFace = 6.0
min_dist_EdgeToFace = 3.8
min_ang_EdgeToFace = 50
max_ang_EdgeToFace = 90

# ==== FaceToFace =============================================================
dist_cut_FaceToFace = 4.5
min_dist_FaceToFace = 3.8
min_ang_FaceToFace = 0
max_ang_FaceToFace = 40

cf = sys.modules[__name__]


def get_cutoff(cutoff_name, args=()):
    """
    Parse the cutoff name from the args or the cf module if not found in args.

    Args:
        args (Namespace): Arguments from the parser.
        cutoff_name (str): Name of the cutoff.

    Returns:
        float: Value of the cutoff.
    """
    if cutoff_name not in cf.__dict__:
        raise ValueError(f"{cutoff_name} is not a valid cutoff name.\n"
                         f"The supported list is:\n"
                         f"{[x for x in dir(cf) if not x.startswith('__')]}")
    elif cutoff_name in args:
        return args[cutoff_name]
    else:
        # get the value from the cf module
        return getattr(cf, cutoff_name)


def parse_cutoffs(args=()):
    """
    Parse the cutoffs from the args or the cf module if not found in args.

    Args:
        args (Namespace): Arguments from the parser.

    Returns:
        dict: Dictionary with the cutoffs.
    """

    cutoffs_raw = {
        # Undirected 1D interactions
        'CloseContacts':
            {'distCut1': get_cutoff('dist_cut_CloseContacts', args)},
        'VdWContact':
            {'distCut1': 0},
        'Hydrophobic':
            {'distCut1': get_cutoff('dist_cut_Hydrophobic', args)},

        # Directed 1D interactions
        'Anionic':
            {'distCut1': get_cutoff('dist_cut_Ionic', args)},
        'Cationic':
            {'distCut1': get_cutoff('dist_cut_Ionic', args)},
        'MetalDonor':
            {'distCut1': get_cutoff('dist_cut_Metalic', args)},
        'MetalAcceptor':
            {'distCut1': get_cutoff('dist_cut_Metalic', args)},

        # Directed 2D1A interactions
        'HBAcceptor':
            {'distCut1': get_cutoff('dist_cut_DA', args),
             'distCut2': get_cutoff('dist_cut_HA', args),
             'minAng1': get_cutoff('min_ang_DHA', args),
             'maxAng1': get_cutoff('max_ang_DHA', args)},
        'HBDonor':
            {'distCut1': get_cutoff('dist_cut_DA', args),
             'distCut2': get_cutoff('dist_cut_HA', args),
             'minAng1': get_cutoff('min_ang_DHA', args),
             'maxAng1': get_cutoff('max_ang_DHA', args)},

        'XBAcceptor':
            {'distCut1': get_cutoff('dist_cut_XA', args),
             'distCut2': get_cutoff('dist_cut_XD', args)},
        'XBDonor':
            {'distCut1': get_cutoff('dist_cut_XA', args),
             'distCut2': get_cutoff('dist_cut_XD', args)},

        # Undirected 2D1A interactions
        'PiStacking':
            {'distCut1': get_cutoff('dist_cut_PiStacking', args),
             'distCut2': get_cutoff('min_dist_PiStacking', args),
             'minAng1': get_cutoff('min_ang_PiStacking', args),
             'maxAng1': get_cutoff('max_ang_PiStacking', args)},
        'FaceToFace':
            {'distCut1': get_cutoff('dist_cut_FaceToFace', args),
             'distCut2': get_cutoff('min_dist_FaceToFace', args),
             'minAng1': get_cutoff('min_ang_FaceToFace', args),
             'maxAng1': get_cutoff('max_ang_FaceToFace', args)},
        'EdgeToFace':
            {'distCut1': get_cutoff('dist_cut_EdgeToFace', args),
             'distCut2': get_cutoff('min_dist_EdgeToFace', args),
             'minAng1': get_cutoff('min_ang_EdgeToFace', args),
             'maxAng1': get_cutoff('max_ang_EdgeToFace', args)},

        # Directed 1D1A interactions
        'PiCation':
            {'distCut1': get_cutoff('dist_cut_PiCation', args),
             'minAng1': get_cutoff('min_ang_PiCation', args),
             'maxAng1': get_cutoff('max_ang_PiCation', args)},

        'CationPi':
            {'distCut1': get_cutoff('dist_cut_PiCation', args),
             'minAng1': get_cutoff('min_ang_PiCation', args),
             'maxAng1': get_cutoff('max_ang_PiCation', args)}}

    return cutoffs_raw


def get_inters_cutoffs(cutoffs_raw):
    """
    Get the interaction names and cutoffs from the raw dictionary

    Args:
        cutoffs_raw (dict): Raw dictionary with the cutoffs.

    Returns:

    """
    inter_names = np.asarray(list(cutoffs_raw.keys()))
    cutoffs_matrix = np.zeros((4, len(inter_names)), dtype=np.float32)

    for i, inter in enumerate(cutoffs_raw):
        if 'distCut1' in cutoffs_raw[inter]:
            cutoffs_matrix[0, i] = cutoffs_raw[inter]['distCut1']

        if 'distCut2' in cutoffs_raw[inter]:
            cutoffs_matrix[1, i] = cutoffs_raw[inter]['distCut2']

        if 'minAng1' in cutoffs_raw[inter]:
            cutoffs_matrix[2, i] = cutoffs_raw[inter]['minAng1']

        if 'maxAng1' in cutoffs_raw[inter]:
            cutoffs_matrix[3, i] = cutoffs_raw[inter]['maxAng1']

    return inter_names, cutoffs_matrix
