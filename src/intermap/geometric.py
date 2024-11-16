# Created by rglez at 11/16/24

import numpy as np
from numba import njit


@njit(parallel=False)
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
    distances = np.empty(n)

    for i in range(n):
        dx = d[i][0] - a[i][0]
        dy = d[i][1] - a[i][1]
        dz = d[i][2] - a[i][2]
        distances[i] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    return distances


@njit(parallel=False)
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
        angle_rad = np.arccos(dot_product / (dh_norm * ah_norm))
        angle_deg[i] = np.rad2deg(angle_rad)  # Convert radians to degrees

    return angle_deg
