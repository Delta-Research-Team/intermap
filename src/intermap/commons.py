# Created by rglez at 12/25/24
"""
Common functions used in the different modules of the package
"""
import logging
import os
import re

import numpy as np
import numpy_indexed as npi
from numba.pycc import CC

cc = CC('commons_aot')

logger = logging.getLogger('InterMapLogger')

smarts = {
    'hydroph': '[c,s,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6&X4&H0]);+0]',
    'cations': '[+{1-},$([N&X3&!$([N&X3]-O)]-C=[N&X3&+])]',
    'anions': '[-{1-},$(O=[C,S,P]-[O&-])]',
    'metal_acc': '[O,#7&!$([n&X3])&!$([N&X3]-*=[!#6])&!$([N&X3]-a)&!$([N&X4]),-{1-};!+{1-}]',
    'metal_don': '[#20,#48,#27,#29,#26,#12,#25,#28,#30]',
    'hb_acc': '[#7&!$([n&X3])&!$([N&X3]-*=[O,N,P,S])&!$([N&X3]-a)&!$([N&v4&+]),O&!$([O&X2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O&-]-N=O),o&+0,F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])]',
    'xb_acc': '[#7,#8,P,S,#34,#52,a;!+{1-}]!#*',
    'hb_don': '[$([O,S;+0]),$([N;v2,v3,v4&+1]),n+0]-[H]',
    'xb_don': '[#6,#7,#14,F,Cl,Br,I]-[Cl,Br,I,#85]',
    'rings5': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:1',
    'rings6': '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:[a&r]:1'}


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

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel("DEBUG")
    logger.addHandler(file_handler)
    return logger


def get_cutoffs_and_inters(to_compute, all_inters, all_cutoffs):
    """
    Get the cutoffs and interactions to compute

    Args:
        to_compute (ndarray): All interactions to compute
        all_inters (ndarray): All interactions
        all_cutoffs (ndarray): All cutoffs

    Returns:
        to_compute_aro (ndarray): Aromatic interactions to compute
        to_compute_others (ndarray): Non-aromatic interactions to compute
        cutoffs_aro (ndarray): Cutoffs for the aromatic interactions
        cutoffs_others (ndarray): Cutoffs for the non-aromatic interactions
    """

    # Parse aromatics
    bit_aro = [y for x in to_compute if re.search(r'Pi|Face', x) for y in
               npi.indices(all_inters, [x])]
    to_compute_aro = np.asarray([all_inters[x] for x in bit_aro])

    # Parse non-aromatics
    bit_others = [y for x in to_compute if not re.search(r'Pi|Face', x) for y
                  in npi.indices(all_inters, [x])]
    to_compute_others = np.asarray([all_inters[x] for x in bit_others],
                                   dtype=str)

    # Get the cutoffs
    cutoffs_aro = all_cutoffs[:, bit_aro]
    cutoffs_others = all_cutoffs[:, bit_others]

    return to_compute_aro, to_compute_others, cutoffs_aro, cutoffs_others


def check_path(path, check_exist=True):
    """
    Check if a path exists and return it if it does

    Args:
        path: path to check
        check_exist: check if the path exists or not

    Returns:
        path: path to the file or directory
    """
    path_exists = os.path.exists(path)
    if check_exist and path_exists:
        return path
    elif (not check_exist) and (not path_exists):
        return path
    elif (not check_exist) and path_exists:
        return path  # todo: check this behaviour
    elif check_exist and (not path_exists):
        raise ValueError(f'\nNo such file or directory: {path}')
    else:
        pass
        raise ValueError(
            f'\nPath already exists and will not be overwritten: {path}')


def check_numeric_in_range(arg_name, value, dtype, minim, maxim):
    """
    Check if a value is of a certain type and within a range

    Args:
        arg_name: name of the argument to check
        value: value to check
        dtype: type of the value
        minim: minimum value
        maxim: maximum value

    Returns:
        value: value of the correct type and within the range
    """
    if not isinstance(value, dtype):
        raise TypeError(f'Param "{arg_name}" must be of type {dtype}')

    if not minim <= value <= maxim:
        raise ValueError(f'Param "{value}" out of [{minim}, {maxim}]')

    return dtype(value)
