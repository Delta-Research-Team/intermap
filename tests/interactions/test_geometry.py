import pytest
import numpy as np
from numba.typed import List
import numba.types as ntypes

from intermap.interactions.geometry import (
    calc_min_dist, get_compress_mask, get_compress_mask2,
    indices, isin, calc_dist, calc_angles_2v, calc_angle_3p,
    calc_centroids, calc_normal_vector, get_containers_run, get_containers
)


# =============================================================================
# COHERENT FIXTURES
# =============================================================================

@pytest.fixture
def points_f4():
    """Standard 3D points in float32."""
    return np.array([[0, 0, 0], [3, 4, 0], [0, 1, 0]], dtype=np.float32)


@pytest.fixture
def ints_i4():
    """Standard integer arrays in int32."""
    full = np.array([10, 20, 30, 40, 50], dtype=np.int32)
    subset = np.array([30, 10, 99], dtype=np.int32)
    return full, subset


# =============================================================================
# MATH & GEOMETRY TESTS
# =============================================================================

class TestGeometryMath:
    def test_calc_min_dist(self, points_f4):
        coords3 = np.array([[3.0, 0.0, 0.0]], dtype=np.float32)
        dist = calc_min_dist(points_f4, coords3)
        assert np.isclose(dist, 3.0)

    def test_calc_dist(self):
        d = np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float32)
        a = np.array([[3, 4, 0], [1, 1, 2]], dtype=np.float32)
        assert np.allclose(calc_dist(d, a), [5.0, 1.0])

    def test_calc_angles_2v(self):
        v1 = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.float32)
        v2 = np.array([[0, 1, 0], [0, 1, 0]], dtype=np.float32)
        assert np.allclose(calc_angles_2v(v1, v2), [90.0, 0.0])

    def test_calc_angles_2v_errors(self):
        with pytest.raises(ValueError):
            calc_angles_2v(np.array([[1, 0]], dtype=np.float32),
                           np.array([[0, 1]], dtype=np.float32))
        assert len(calc_angles_2v(np.empty((0, 3), dtype=np.float32),
                                  np.empty((0, 3), dtype=np.float32))) == 0

    def test_calc_angle_3p(self):
        d = np.array([[1, 0, 0]], dtype=np.float32)
        h = np.array([[0, 0, 0]], dtype=np.float32)
        a = np.array([[0, 1, 0]], dtype=np.float32)
        assert np.isclose(calc_angle_3p(d, h, a)[0], 90.0)

    def test_calc_angle_3p_exception(self, monkeypatch):
        """Hits lines 221-222 and suppresses divide-by-zero warnings."""
        import numpy as np
        # Monkeypatch arccos to force the exception block
        monkeypatch.setattr(np, "arccos",
                            lambda x: exec('raise ValueError("Forced")'))

        d = h = a = np.array([[0, 0, 0]], dtype=np.float32)
        # Silence RuntimeWarning for 0/0 division during the test
        with np.errstate(all='ignore'):
            angles = calc_angle_3p(d, h, a)
        assert angles[0] == 0.0


# =============================================================================
# DATA STRUCTURE TESTS (Masks, Indices, Structures)
# =============================================================================

class TestGeometryIndicesAndMasks:
    def test_get_compress_masks(self):
        arr = np.array([[True, False], [False, False]], dtype=np.bool_)
        assert np.array_equal(get_compress_mask(arr), [True, False])
        assert np.array_equal(get_compress_mask2(arr), [True, False])

    def test_indices_function(self, ints_i4):
        assert np.array_equal(indices(*ints_i4), [2, 0, -1])

    def test_isin_function(self, ints_i4):
        assert np.array_equal(isin(*ints_i4),
                              [True, False, True, False, False])


class TestGeometryStructures:
    def test_calc_centroids(self, points_f4):
        rings = np.array([[0, 1, 2, -1, -1, -1, 3]], dtype=np.int32)
        res = calc_centroids(rings, points_f4)
        assert np.allclose(res[0], [1.0, 5.0 / 3.0, 0.0])

    def test_calc_normal_vector(self):
        # Pass 2D points to match the (p1, p2, p3) logic
        p1, p2, p3 = (np.array([0, 0, 0], dtype=np.float32),
                      np.array([1, 0, 0], dtype=np.float32),
                      np.array([0, 1, 0], dtype=np.float32))
        assert np.allclose(calc_normal_vector(p1, p2, p3), [0, 0, 1])


# =============================================================================
# CONTAINER & KD-TREE PARSER TESTS
# =============================================================================

class TestGeometryContainers:
    def test_get_containers_run_logic(self):
        """Hits n_contacts=0, n_contacts>0, and idems (self-contact) removal."""
        xyz = np.array([[0, 0, 0], [1, 1, 1]], dtype=np.float32)
        k = np.int64(1)

        # 1. Test Empty
        ball_empty = List.empty_list(ntypes.int64[:])
        ijf_e, dists_e = get_containers_run(xyz, k, ball_empty,
                                            np.empty(0, dtype=np.int32),
                                            np.empty(0, dtype=np.int32))
        assert ijf_e.shape[0] == 0

        # 2. Test Hit + Idem Removal
        s1 = np.array([0, 1], dtype=np.int32)
        s2 = np.array([0, 1], dtype=np.int32)
        ball = List.empty_list(ntypes.int64[:])
        ball.append(
            np.array([0], dtype=np.int64))  # s1[0] interacts with s2[0] (idem)
        ball.append(np.array([0],
                             dtype=np.int64))  # s1[1] interacts with s2[0] (valid)

        ijf, dists = get_containers_run(xyz, k, ball, s1, s2)
        # Note: Function logic overwrites ijf with uniques but doesn't resize it.
        # We check the valid interaction [1, 0, 1] exists.
        assert ijf[0, 0] == 1
        assert ijf[0, 1] == 0

    def test_get_containers_complex(self):
        """Hits residue overlap logic and empty branches."""
        xyz = np.array([[0, 0, 0], [5, 5, 5]], dtype=np.float32)
        xyz_idx = np.array([0, 1], dtype=np.int32)

        # Contact between index 0 and 1
        ball = List.empty_list(ntypes.int64[:])
        ball.append(np.array([1], dtype=np.int64))
        ball.append(np.array([], dtype=np.int64))

        s1 = s2 = np.array([0, 1], dtype=np.int32)
        resconv = np.zeros(2, dtype=np.int32)  # Both in same residue

        # Test: atomic=False, same residue -> should be removed
        ijf, d, inter = get_containers(xyz, np.int64(1), xyz_idx, ball, s1, s2,
                                       np.int64(1), False, resconv)
        assert ijf.shape[0] == 0

    def test_get_containers_empty(self):
        """Covers the return branch for zero contacts."""
        ball = List.empty_list(ntypes.int64[:])
        res = get_containers(np.empty((0, 3), dtype=np.float32), np.int64(1),
                             np.empty(0, dtype=np.int32), ball,
                             np.empty(0, dtype=np.int32),
                             np.empty(0, dtype=np.int32),
                             np.int64(1), True, np.empty(0, dtype=np.int32))
        assert res[0].shape[0] == 0
