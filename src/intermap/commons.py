# Created by rglez at 12/25/24
"""
Common functions used in the different modules of the package
"""
import logging
import os
import re

import numpy as np
import numpy_indexed as npi
from numba import njit, prange

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


@njit(parallel=False, cache=True)
def get_compress_mask(array):
    """
    Compress the array to remove empty

    Args:
        array (ndarray): Array to compress

    Returns:
        array (ndarray): Compressed array
    """
    n = array.shape[0]
    mask = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        if array[i].any():
            mask[i] = True
    return mask


@njit(parallel=False, cache=True)
def indices(full, subset):
    """
    Find the indices of the subset elements in the full array.

    Args:
        full (ndarray): Full array.
        subset (ndarray): Subset array.

    Returns:
        indices (ndarray): Indices of the subset elements in the full array.
    """
    indices = np.full(len(subset), -1, dtype=np.int32)
    for i in range(len(subset)):
        for j in range(len(full)):
            if full[j] == subset[i]:
                indices[i] = j
                break
    return indices


@njit(parallel=False, cache=True)
def isin(full, subset):
    """
    Check if the elements of the subset are in the full array.

    Args:
        full (ndarray): Full array.
        subset (ndarray): Subset array.

    Returns:
        result (ndarray): Boolean array indicating if the elements of the
                          subset are in the full array.
    """
    n = len(full)
    result = np.full(n, False)
    set_b = set(subset)
    for i in prange(n):
        if full[i] in set_b:
            result[i] = True
    return result


@njit(parallel=False, cache=True)
def calc_dist(d, a):
    """
    Computes the Euclidean distance between two atoms in a molecule

    Args:
        d (ndarray): Coordinates of the first atom (n, 3).
        a (ndarray): Coordinates of the second atom (n, 3).

    Returns:
        float: the Euclidean distance between the two atoms
    """
    n = d.shape[0]
    distances = np.empty(n, dtype=np.float32)

    for i in prange(n):
        dx = d[i][0] - a[i][0]
        dy = d[i][1] - a[i][1]
        dz = d[i][2] - a[i][2]
        distances[i] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return distances


@njit(parallel=False, cache=True)
def calc_min_dist(coords1, coords2):
    """
    Get the minimumm distance between two sets of coordinates

    Args:
        coords1: coordinates of the first residue
        coords2: coordinates of the second residue

    Returns:
        The minimum distance between two sets of coordinates
    """
    # Constants
    n1 = coords1.shape[0]
    n2 = coords2.shape[0]

    # Find minimum distance using square values to save time
    min_dist_squared = np.inf
    for i in range(n1):
        for j in range(n2):
            dist_squared = (
                    (coords1[i][0] - coords2[j][0]) ** 2 +
                    (coords1[i][1] - coords2[j][1]) ** 2 +
                    (coords1[i][2] - coords2[j][2]) ** 2)
            if dist_squared < min_dist_squared:
                min_dist_squared = dist_squared
    return np.sqrt(min_dist_squared)


@njit(parallel=False, cache=True)
def calc_angles_2v(vectors1, vectors2):
    """
    Computes the angles between two arrays of 3D vectors using arctan2.

    Args:
        vectors1 (np.ndarray): Array of shape (N, 3) with the first set of vectors.
        vectors2 (np.ndarray): Array of shape (N, 3) with the second set of vectors.

    Returns:
        np.ndarray: Array of shape (N,) containing the angles in degrees.
    """

    # Input checks and validations
    have_same_shape = vectors1.shape == vectors2.shape
    are_2d_arrays = vectors1.ndim == 2 and vectors2.ndim == 2
    have_3_elements = vectors1.shape[1] == 3
    sms = "Input arrays must be 2D arrays of shape (N, 3) and have the same shape."

    if not (have_same_shape and are_2d_arrays and have_3_elements):
        raise ValueError(sms)
    if vectors1.size == 0:
        return None

    num_vectors = vectors1.shape[0]
    angles_degrees = np.empty(num_vectors)

    for i in range(num_vectors):
        # Compute dot product
        v1 = vectors1[i]
        v2 = vectors2[i]
        dot_product = np.sum(v1 * v2)

        # Compute cross product
        cross_product = np.empty(3)
        cross_product[0] = v1[1] * v2[2] - v1[2] * v2[1]
        cross_product[1] = v1[2] * v2[0] - v1[0] * v2[2]
        cross_product[2] = v1[0] * v2[1] - v1[1] * v2[0]

        # Compute angle
        cross_product_magnitude = np.sqrt(np.sum(cross_product ** 2))
        angle_radians = np.arctan2(cross_product_magnitude, dot_product)
        angles_degrees[i] = np.degrees(angle_radians)

    return angles_degrees


@njit(parallel=False, cache=True)
def calc_angle(d, h, a):
    """
    Computes the angles between sets of three atoms.

    Args:
        d (ndarray): Coordinates of the donor atoms (n, 3).
        h (ndarray): Coordinates of the hydrogen atoms (n, 3).
        a (ndarray): Coordinates of the acceptor atoms (n, 3).

    Returns:
        angle_deg: Angles in degrees for each set of atoms (n,).
    """
    n = d.shape[0]  # Number of triplets
    angle_deg = np.empty(n)  # Array to hold the result angles

    for i in range(n):
        # Compute vectors
        dh = d[i] - h[i]
        ah = a[i] - h[i]

        # Compute dot product and norms
        dot_product = dh[0] * ah[0] + dh[1] * ah[1] + dh[2] * ah[2]
        dh_norm = np.sqrt(dh[0] ** 2 + dh[1] ** 2 + dh[2] ** 2)
        ah_norm = np.sqrt(ah[0] ** 2 + ah[1] ** 2 + ah[2] ** 2)

        # Compute angle
        try:
            angle_rad = np.arccos(dot_product / (dh_norm * ah_norm))
            angle_deg[i] = np.rad2deg(angle_rad)  # Convert radians to degrees
        except:
            angle_deg[i] = 0.0

    return angle_deg


@njit(parallel=False, cache=True)
def calc_centroids(rings, xyz):
    """
    Calculate the centroid of each ring

    Args:
        rings (list): List of rings
        xyz (np.ndarray): Coordinates of the atoms

    Returns:
        centroids (np.ndarray): Centroids of the rings
    """
    centroids = np.zeros((len(rings), 3), dtype=float)
    for i, ring in enumerate(rings):
        atoms = ring[:ring[-1]]
        centroids[i] = xyz[atoms].sum(axis=0) / len(atoms)
    return centroids.astype(np.float32)


@njit(parallel=False, cache=True)
def calc_normal_vector(p1, p2, p3):
    """
    Calculate the normal vector of a plane defined by three points

    Args:
        p1 (np.ndarray): First point
        p2 (np.ndarray): Second point
        p3 (np.ndarray): Third point

    Returns:
        np.ndarray: Normal vector of the plane defined by the three points
    """

    # Calculate vectors from points
    v1 = p2 - p1
    v2 = p3 - p1

    # Calculate the normal vector
    normal = np.cross(v1, v2)
    norm = np.linalg.norm(normal)

    return (normal / norm).astype(np.float32)


@njit(parallel=False, cache=True)
def get_containers(xyz, k, ext_idx, ball_1, s1_indices, s2_indices,
                   to_compute):
    """
    Get the containers for the interactions, the distances and the indices

    Args:
        xyz (ndarray): Coordinates of the atoms in the system
        k (int): Index of the frame in the trajectory
        ext_idx (ndarray): External (real) indices for the atoms in the system
        ball_1 (list): List of lists with the indices of the atoms in the second selection
        s1_indices (ndarray): Indices of the atoms in the first selection
        s2_indices (ndarray): Indices of the atoms in the second selection
        to_compute (ndarray): Interactions to compute

    Returns:
        ijf (ndarray): Indices of the atoms in the first and second selections
        dists (ndarray): Distances between the atoms in the first and second selections
        interactions (ndarray): Container for the interactions
    """

    # Find the contacts
    n_contacts = sum([len(x) for x in ball_1])
    ijf = np.zeros((n_contacts, 3), dtype=np.int32)
    counter = 0
    for i, x in enumerate(ball_1):
        X1 = s1_indices[i]
        for j in x:
            X2 = s2_indices[j]
            ijf[counter][0] = X1
            ijf[counter][1] = X2
            ijf[counter][2] = k
            counter += 1

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = ext_idx[ijf[:, 0]] == ext_idx[ijf[:, 1]]
    if idems.any():
        ijf = ijf[~idems]

    # Compute distances
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    n_types = to_compute.size
    interactions = np.zeros((ijf.shape[0], n_types), dtype=np.bool_)
    return ijf, dists, interactions
