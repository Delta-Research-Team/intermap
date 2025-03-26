# Created by rglez at 12/25/24
"""
Common functions used in the different modules of the package
"""
import logging
import os

import numpy as np
from numba import njit, prange
from numba.typed import Dict, List
from numba_kdtree import KDTree as nckd

logger = logging.getLogger('InterMapLogger')


def trajiter(universe, chunk_frames, sel_idx):
    for chunk in chunk_frames:
        xyz_chunk = np.empty((chunk.size, sel_idx.size, 3),
                             dtype=np.float32)
        for i, frame in enumerate(chunk):
            xyz_chunk[i] = universe.trajectory[frame].positions[
                sel_idx]
        yield xyz_chunk


def parse_last_param(last, traj_len):
    """
    Parse the last parameter

    Args:
        last (int): Last frame to process
        traj_len (int): Length of the trajectory

    Returns:
        last (int): Last frame to process
    """
    if last == -1:
        last = traj_len
    elif last > traj_len:
        logger.warning(
            f"The declared end frame '{last}' is larger than the trajectory "
            f"length '{traj_len}'. The end frame will be set to the last frame.")
        last = traj_len
    return last


def split_in_chunks(array, chunk_size):
    """
    Split an array in chunks of the specified size

    Args:
        array (ndarray): Array to split
        chunk_size (int): Size of the chunks

    Returns:
        chunks (generator): Generator with the chunks
    """
    for i in range(0, array.shape[0], chunk_size):
        yield array[i:i + chunk_size]


def get_coordinates(u, chunk, sel_idx):
    xyz_chunk = np.empty((chunk.size, sel_idx.size, 3), dtype=np.float32)
    for i, frame in enumerate(chunk):
        xyz_chunk[i] = u.trajectory[frame].positions[sel_idx].copy()
    return xyz_chunk


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


def get_trees(xyz_chunk, s2_indices):
    """

    Args:
        xyz_chunk:
        s2_indices:

    Returns:

    """
    trees = List()
    for x in xyz_chunk:
        trees.append(nckd(x[s2_indices]))
    return trees


@njit(parallel=True, cache=True)
def get_trees2(xyz_chunk, s2_idx):
    """

    Args:
        xyz_chunk:
        s2_idx:

    Returns:

    """
    N = len(xyz_chunk)
    trees = Dict()
    value = nckd(xyz_chunk[0][s2_idx])
    for i in range(N):
        trees[i] = value

    for i in prange(N):
        trees[i] = nckd(xyz_chunk[i][s2_idx])
    return trees


@njit(parallel=False, cache=True)
def get_ball(xyz, s1_indices, tree, dist_cut):
    """
    Args:
        xyz:
        s1_indices:
        tree:
        dist_cut:

    Returns:

    """
    ball_1 = tree.query_radius(xyz[s1_indices], dist_cut)
    return ball_1
