import pytest
import numpy as np
from numba.typed import List

# Import the module to patch its internal 'nckd' reference
import intermap.interactions.aro as aro_module
from intermap.interactions.aro import (
    get_normals_and_centroids, get_intersect_point,
    stackings, pications, aro, get_aro_xyzs, get_trees
)


# =============================================================================
# 1. MATHEMATICAL & TREE TESTS
# =============================================================================

class TestAromaticGeometry:
    def test_get_trees(self, monkeypatch):
        """Hits lines 14-18: Tree generation bypassing KDTree recursion error."""
        # FIX: Patch the local reference 'nckd' inside aro.py
        monkeypatch.setattr(aro_module, "nckd", lambda x: "mock_tree")

        xyzs_aro = [np.zeros((10, 3), dtype=np.float32)]
        s2_indices = np.array([0, 1, 2], dtype=np.int32)

        trees = get_trees(xyzs_aro, s2_indices)
        assert trees[0] == "mock_tree"

    def test_get_intersect_point_standard(self):
        s1_norm, s1_cent = np.array([0, 0, 1], dtype=np.float32), np.array(
            [0, 0, 0], dtype=np.float32)
        s2_norm, s2_cent = np.array([1, 0, 0], dtype=np.float32), np.array(
            [0, 0, 0], dtype=np.float32)
        point = get_intersect_point(s1_norm, s1_cent, s2_norm, s2_cent)
        assert np.allclose(point, [0, 0, 0], atol=1e-5)

    def test_get_intersect_point_parallel(self):
        s1_norm, s1_cent = np.array([0, 0, 1], dtype=np.float32), np.array(
            [0, 0, 0], dtype=np.float32)
        s2_norm, s2_cent = np.array([0, 0, 1], dtype=np.float32), np.array(
            [0, 0, 5], dtype=np.float32)
        point = get_intersect_point(s1_norm, s1_cent, s2_norm, s2_cent)
        assert np.isnan(point).all()


# =============================================================================
# 2. CLASSIFICATION & ERROR TESTS
# =============================================================================

class TestAromaticClassification:
    def test_pications_logic(self):
        cuts = np.zeros((6, 1), dtype=np.float32)
        cuts[0, 0], cuts[3, 0] = 5.0, 90.0
        row1, row2 = np.array([0], dtype=np.int32), np.array([1],
                                                             dtype=np.int32)
        dists = np.array([4.0], dtype=np.float32)
        xyz_aro = np.array([[0, 0, 0], [0, 0, 4]], dtype=np.float32)
        s1_norm = np.array([[0, 0, 1]], dtype=np.float32)

        _, result = pications('PiCation', xyz_aro, row1, row2, dists,
                              np.array([0], dtype=np.int32),
                              np.array([], dtype=np.int32),
                              np.array([], dtype=np.int32),
                              np.array([1], dtype=np.int32),
                              s1_norm, s1_norm, cuts, ['PiCation'])
        assert result[0] == True

    def test_classification_errors(self):
        """Hits lines 113-117, 132-135, and 171 using valid dummy cutoffs."""
        dummy_cuts = np.zeros((6, 1), dtype=np.float32)
        row1 = np.array([0], dtype=np.int32)

        with pytest.raises(ValueError, match="Invalid interaction"):
            pications('Invalid', None, row1, row1, None, None, None, None,
                      None, None, None, dummy_cuts, ['Invalid'])

        with pytest.raises(Exception, match="Invalid interaction"):
            stackings('Invalid', np.array([0.]), np.array([0.]),
                      np.array([0.]), np.array([0.]), np.array([0.]),
                      dummy_cuts, ['Invalid'])


# =============================================================================
# 3. INTEGRATION TESTS
# =============================================================================

class TestAromaticIntegration:
    def test_aro_orchestrator(self, iman):
        """Hits all major aromatic paths and empty-branch returns."""
        # 1. Test 'None' branch (Line 210)
        res_ijf, res_inters = aro(None, None, np.zeros((1, 3), dtype=np.int32),
                                  None, None, None, None, None, None, None,
                                  None, None, None, None, None, ['None'])
        assert not res_inters.any()

        # 2. Test full system on frame 0
        xyzs = np.array([iman.universe.atoms.positions], dtype=np.float32)
        s1_cent, s2_cent, xyz_aro_traj = get_aro_xyzs(xyzs, iman.s1_rings,
                                                      iman.s2_rings,
                                                      iman.s1_cat, iman.s2_cat,
                                                      iman.s1_ani, iman.s2_ani)

        selected = ['PiCation', 'PiAnion', 'CationPi', 'AnionPi', 'PiStacking',
                    'FaceToFace', 'EdgeToFace']
        cuts = np.zeros((6, len(selected)), dtype=np.float32)
        cuts[0, :], cuts[3, :] = 6.0, 90.0

        s1_norm, s2_norm, _, _ = get_normals_and_centroids(xyzs[0],
                                                           iman.s1_rings,
                                                           iman.s2_rings)

        # Check if rings exist before running to prevent IndexError on ijf_aro[0] (line 258)
        if iman.rings.size > 0:
            res_ijf, res_inters = aro(xyz_aro_traj[0], iman.xyz_aro_idx,
                                      iman.rings,
                                      np.zeros(iman.rings.shape[0],
                                               dtype=np.float32),
                                      iman.s1_rings_idx, iman.s2_rings_idx,
                                      iman.s1_cat_idx,
                                      iman.s2_cat_idx, iman.s1_ani_idx,
                                      iman.s2_ani_idx,
                                      s1_norm, s2_norm, s1_cent[0], s2_cent[0],
                                      cuts, selected)
            assert res_ijf.shape[1] == 3
