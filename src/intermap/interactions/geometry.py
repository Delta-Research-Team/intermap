# Created by rglez at 2/8/25
"""
Common functions used in the different modules of the package
"""
import numpy as np
from numba import njit, prange


@njit("f4(f4[:, :], f4[:, :])", parallel=False, cache=True)
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
    return np.float32(np.sqrt(min_dist_squared))


@njit("b1[:](b1[:, :])", parallel=False, cache=True)
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
def get_compress_mask2(array):
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


@njit("i4[:](i4[:], i4[:])", parallel=False, cache=True)
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


@njit("b1[:](i4[:], i4[:])", parallel=False, cache=True)
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


@njit("f4[:](f4[:, :], f4[:, :])", parallel=False, cache=True)
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


@njit("f4[:](f4[:, :], f4[:, :])", parallel=False, cache=True)
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
        return np.zeros(0, dtype=np.float32)

    num_vectors = vectors1.shape[0]
    angles_degrees = np.empty(num_vectors, dtype=np.float32)

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


@njit("f4[:](f4[:, :], f4[:, :], f4[:, :])", parallel=False, cache=True)
def calc_angle_3p(d, h, a):
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
    angle_deg = np.empty(n,
                         dtype=np.float32)  # Array to hold the result angles

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


@njit("f4[:, :](i4[:, :], f4[:, :])", parallel=False, cache=True)
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


@njit("f4[:, :](f4[:, :], f4[:, :], f4[:, :])", parallel=False, cache=True)
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


@njit("Tuple((i4[:, :], f4[:]))"
      "(f4[:, :], i8, ListType(i8[:]), i4[:], i4[:])",
      cache=True)
def get_containers_run(xyz, k, ball_aro, s1_idx, s2_idx, ):
    # Find the contacts
    n_contacts = sum([len(x) for x in ball_aro])
    ijf = np.zeros((n_contacts, 3), dtype=np.int32)

    if n_contacts == 0:
        return ijf, np.zeros(0, dtype=np.float32)

    counter = 0
    for i, x in enumerate(ball_aro):
        X1 = s1_idx[i]
        for j in x:
            X2 = s2_idx[j]
            ijf[counter][0] = X1
            ijf[counter][1] = X2
            ijf[counter][2] = k
            counter += 1

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = ijf[:counter, 0] == ijf[:counter, 1]
    uniques = ijf[:counter][~idems].copy()
    n_uniques = uniques.shape[0]
    ijf[:n_uniques] = uniques

    # Compute distances
    row1 = ijf[:n_uniques, 0]
    row2 = ijf[:n_uniques, 1]
    dists = calc_dist(xyz[row1], xyz[row2])

    # Create the container for interaction types
    return ijf, dists


@njit("f4[:, :](f4[:, :], f4[:, :], f4[:, :])", parallel=False, cache=True)
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


@njit("Tuple((i4[:, :], f4[:], b1[:, :]))"
      "(f4[:, :], i8, i4[:], ListType(i8[:]), i4[:], i4[:], i8, b1, i4[:])",
      cache=True)
def get_containers(xyz_aro, k, xyz_aro_idx, ball_aro, s1_aro_idx, s2_aro_idx,
                   n_types, atomic, resconv):
    # Find the contacts
    n_contacts = sum([len(x) for x in ball_aro])
    if n_contacts == 0:
        return (
            np.zeros((0, 3), dtype=np.int32), np.zeros(0, dtype=np.float32),
            np.zeros((0, n_types), dtype=np.bool_))

    ijf = np.zeros((n_contacts, 3), dtype=np.int32)
    counter = 0
    for i, x in enumerate(ball_aro):
        X1 = s1_aro_idx[i]
        for j in x:
            X2 = s2_aro_idx[j]
            ijf[counter][0] = X1
            ijf[counter][1] = X2
            ijf[counter][2] = k
            counter += 1

    # Remove self-residues if not atomic resolution requested
    if not atomic:
        r1 = resconv[xyz_aro_idx[ijf[:, 0]]]
        r2 = resconv[xyz_aro_idx[ijf[:, 1]]]
        idems = r1 == r2
        ijf = ijf[~idems]

    # Remove idems (self-contacts appearing if both selections overlap)
    idems = xyz_aro_idx[ijf[:, 0]] == xyz_aro_idx[ijf[:, 1]]
    if idems.any():
        ijf = ijf[~idems]

    # Compute distances
    row1 = ijf[:, 0]
    row2 = ijf[:, 1]
    dists = calc_dist(xyz_aro[row1], xyz_aro[row2])

    # Create the container for interaction types
    interactions = np.zeros((ijf.shape[0], n_types), dtype=np.bool_)
    return ijf, dists, interactions
