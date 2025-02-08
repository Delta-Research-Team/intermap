# Created by gonzalezroy at 2/6/25
"""
Test the geometric functions to compute distances and angles
"""
import mdtraj as md
import numpy as np
import pytest

import intermap.njitted as aot
from tests.conftest import mdtrajectory


def test_calc_normal_vector_valid_input():
    p1 = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float32)
    p2 = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    p3 = np.array([[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32)
    result = aot.calc_normal_vector(p1, p2, p3)
    expected = np.array([[0.0, 0.0, 0.70710677], [0.0, 0.0, 0.70710677]],
                        dtype=np.float32)
    assert np.allclose(result, expected)


@pytest.mark.xfail
def test_calc_centroids_valid_input():
    rings = np.array([[0, 1, 2, 3], [4, 5, 6, 7]], dtype=np.int32)
    xyz = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0]
    ], dtype=np.float32)
    result = aot.calc_centroids(rings, xyz)
    expected = np.array([
        [0.5, 0.5, 0.0],
        [0.5, 0.5, 1.0]
    ], dtype=np.float32)
    assert np.allclose(result, expected)


def test_calc_angle_valid_input():
    d = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32)
    h = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], dtype=np.float32)
    a = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float32)
    result = aot.calc_angle(d, h, a)
    expected = np.array([90.0, 90.0], dtype=np.float32)
    assert np.allclose(result, expected)


def test_calc_angles_2v_valid_input():
    vectors1 = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.float32)
    vectors2 = np.array([[0, 1, 0], [1, 0, 0]], dtype=np.float32)
    result = aot.calc_angles_2v(vectors1, vectors2)
    expected = np.array([90.0, 90.0], dtype=np.float32)
    assert np.allclose(result, expected)


def test_calc_dist_valid_input():
    d = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float32)
    a = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]], dtype=np.float32)
    dist = np.linalg.norm(d - a, axis=1)
    result = aot.calc_dist(d, a)
    assert np.allclose(result, dist)


def test_isin_valid_input():
    full = np.array([10, 20, 30, 40, 50], dtype=np.int32)
    subset = np.array([20, 40], dtype=np.int32)
    result = aot.isin(full, subset)
    expected = np.array([False, True, False, True, False], dtype=np.bool_)
    assert np.array_equal(result, expected)


def test_indices_valid_input():
    full = np.array([10, 20, 30, 40, 50], dtype=np.int32)
    subset = np.array([20, 40], dtype=np.int32)
    result = aot.indices(full, subset)
    expected = np.array([1, 3], dtype=np.int32)
    assert np.array_equal(result, expected)


def test_get_compress_mask_valid_input():
    array = np.array([[1, 0], [0, 0], [1, 1]], dtype=bool)
    result = aot.get_compress_mask(array)
    expected = np.array([True, False, True], dtype=bool)
    assert np.array_equal(result, expected)


def get_compress_mask_empty_input():
    array = np.array([], dtype=bool).reshape(0, 2)
    result = aot.get_compress_mask(array)
    expected = np.array([], dtype=bool)
    assert np.array_equal(result, expected)


def get_compress_mask_all_empty_rows():
    array = np.array([[0, 0], [0, 0], [0, 0]], dtype=bool)
    result = aot.get_compress_mask(array)
    expected = np.array([False, False, False], dtype=bool)
    assert np.array_equal(result, expected)


def get_compress_mask_all_non_empty_rows():
    array = np.array([[1, 0], [0, 1], [1, 1]], dtype=bool)
    result = aot.get_compress_mask(array)
    expected = np.array([True, True, True], dtype=bool)
    assert np.array_equal(result, expected)


def get_compress_mask_mixed_rows():
    array = np.array([[1, 0], [0, 0], [0, 1]], dtype=bool)
    result = aot.get_compress_mask(array)
    expected = np.array([True, False, True], dtype=bool)
    assert np.array_equal(result, expected)


def test_min_dist_valid_input():
    coords1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]], dtype=np.float32)
    coords2 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]], dtype=np.float32)
    result = aot.calc_min_dist(coords1, coords2)
    assert result == np.sqrt(1.0)


def test_min_dist_identical_points():
    coords1 = np.array([[1.0, 1.0, 1.0]], dtype=np.float32)
    coords2 = np.array([[1.0, 1.0, 1.0]], dtype=np.float32)
    result = aot.calc_min_dist(coords1, coords2)
    assert result == 0.0


def test_distance(mdtrajectory):
    """
    Test the distance function
    """

    trj = mdtrajectory
    n_atoms = trj.n_atoms

    # get random indices for the tests
    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)

    # compute distances with mdtraj
    dists_mdtraj = md.compute_distances(trj, zip(idx1, idx2), periodic=False)

    # compute distances with imap
    xyz1 = trj.xyz[:, idx1][0]
    xyz2 = trj.xyz[:, idx2][0]
    dists_imap = aot.calc_dist(xyz1, xyz2)

    assert np.allclose(dists_mdtraj, dists_imap, atol=1e-5)


def test_angles_3p(mdtrajectory):
    trj = mdtrajectory
    n_atoms = trj.n_atoms

    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)
    idx3 = np.random.randint(0, n_atoms, 500)

    abc = np.stack([idx1, idx2, idx3], axis=1)
    no_equals = (abc[:, 0] != abc[:, 1]) & (abc[:, 1] != abc[:, 2]) & (
            abc[:, 0] != abc[:, 2])
    trio = abc[no_equals]

    # compute angles with mdtraj
    angles_mdtraj_radians = md.compute_angles(trj, trio, periodic=False)
    angles_mdtraj = np.degrees(angles_mdtraj_radians)

    # compute angles with imap
    xyz1 = trj.xyz[:, trio[:, 0]][0]
    xyz2 = trj.xyz[:, trio[:, 1]][0]
    xyz3 = trj.xyz[:, trio[:, 2]][0]
    angles_imap = aot.calc_angle(xyz1, xyz2, xyz3)

    assert np.allclose(angles_mdtraj, angles_imap, atol=1e-5)


def test_angles_2v(mdtrajectory):
    trj = mdtrajectory
    n_atoms = trj.n_atoms

    idx1 = np.random.randint(0, n_atoms, 500)
    idx2 = np.random.randint(0, n_atoms, 500)
    idx3 = np.random.randint(0, n_atoms, 500)

    abc = np.stack([idx1, idx2, idx3], axis=1)
    no_equals = (abc[:, 0] != abc[:, 1]) & (abc[:, 1] != abc[:, 2]) & (
            abc[:, 0] != abc[:, 2])
    trio = abc[no_equals]

    # compute angles with mdtraj
    angles_mdtraj_radians = md.compute_angles(trj, trio, periodic=False)
    angles_mdtraj = np.degrees(angles_mdtraj_radians)[0]

    # compute angles with imap 2v
    xyz1 = trj.xyz[:, trio[:, 0]][0]
    xyz2 = trj.xyz[:, trio[:, 1]][0]
    xyz3 = trj.xyz[:, trio[:, 2]][0]

    v1 = xyz1 - xyz2
    v2 = xyz3 - xyz2
    angles_imap_2v = aot.calc_angles_2v(v1, v2)
    angles_imap_3p = aot.calc_angle(xyz1, xyz2, xyz3)

    assert np.allclose(angles_mdtraj, angles_imap_2v, atol=1e-3)
    assert np.allclose(angles_mdtraj, angles_imap_3p, atol=1e-3)
    assert np.allclose(angles_imap_2v, angles_imap_3p, atol=1e-3)
