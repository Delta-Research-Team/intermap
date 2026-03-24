import pytest
import numpy as np
from numba.typed import List
import intermap.interactions.others as others_module
from intermap.interactions.others import (
    unswap_frame, detect_vdw, detect_1d, detect_hbonds, others, containers,
    get_ball
)


# =============================================================================
# DUMMY MATH FIXTURES
# =============================================================================

@pytest.fixture
def dummy_coords():
    return np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [10.0, 10.0, 10.0]
    ], dtype=np.float32)


# =============================================================================
# UNIT TESTS (Logic & Math)
# =============================================================================

class TestOthersMathUnit:
    def test_unswap_frame_logic(self):
        ijf = np.array([[0, 1, 0], [1, 0, 0], [2, 3, 0]], dtype=np.int32)
        assert unswap_frame(ijf).sum() == 2

    def test_detect_vdw_math(self):
        dists = np.array([2.0, 2.0], dtype=np.float32)
        row1, row2 = np.array([0, 2], dtype=np.int32), np.array([1, 3],
                                                                dtype=np.int32)
        radii = np.array([1.05, 1.05, 0.7, 0.7], dtype=np.float32)
        _, mask = detect_vdw(dists, row1, row2, radii, ['VdWContact'])
        assert mask[0] == True and mask[1] == False

    def test_detect_1d_logic(self):
        dists = np.array([3.0, 6.0], dtype=np.float32)
        row1, row2 = np.array([0, 1], dtype=np.int32), np.array([1, 0],
                                                                dtype=np.int32)
        type_arr = np.array([0, 1], dtype=np.int32)
        _, m1 = detect_1d('CloseContact', dists, row1, type_arr, row2,
                          type_arr, np.array([[5.0]], dtype=np.float32),
                          ['CloseContact'])
        assert m1[0] == True

    def test_detect_hbonds_donor_math(self, dummy_coords):
        row1, row2 = np.array([1], dtype=np.int32), np.array([2],
                                                             dtype=np.int32)
        _, mask = detect_hbonds("HBDonor", row1, row1, row2, row2,
                                np.array([2.0], dtype=np.float32),
                                dummy_coords, np.zeros(4, dtype=np.int32), 2.5,
                                3.5, 120.0, 180.0, ["HBDonor"])
        assert mask[0] == True

    def test_detect_xbonds_acceptor_math(self, dummy_coords):
        row1, row2 = np.array([2], dtype=np.int32), np.array([1],
                                                             dtype=np.int32)
        _, mask = detect_hbonds("XBAcceptor", row1, row1, row2, row2,
                                np.array([2.0], dtype=np.float32),
                                dummy_coords, np.zeros(4, dtype=np.int32), 2.5,
                                3.5, 120.0, 180.0, ["XBAcceptor"])
        assert mask[0] == True

    def test_detect_hbonds_empty_and_invalid(self, dummy_coords):
        row1, row2 = np.array([1], dtype=np.int32), np.array([2],
                                                             dtype=np.int32)
        # Miss distance cutoff
        _, mask = detect_hbonds("HBDonor", row1, row1, row2, row2,
                                np.array([10.0], dtype=np.float32),
                                dummy_coords, np.zeros(4, dtype=np.int32), 2.5,
                                3.5, 120.0, 180.0, ["HBDonor"])
        assert not mask.any()
        # Invalid Name
        with pytest.raises(ValueError, match="Invalid interaction name"):
            detect_hbonds("Invalid", row1, row1, row2, row2,
                          np.array([2.0], dtype=np.float32),
                          dummy_coords, np.zeros(4, dtype=np.int32), 2.5, 3.5,
                          120.0, 180.0, ["Invalid"])

    def test_get_ball_mock(self, monkeypatch):
        class MockTree:
            def query_radius(self, x, r): return [
                np.array([1], dtype=np.int64)]

        monkeypatch.setattr(others_module, "nckd", lambda x: MockTree())
        assert len(
            get_ball(np.zeros((2, 3), dtype=np.float32), [0], [1], 5.0)) == 1


# =============================================================================
# ORCHESTRATOR TESTS
# =============================================================================

class TestOthersOrchestrator:
    def test_hb_xb_partial_branches(self):
        xyz = np.zeros((2, 3), dtype=np.float32)
        ball = [np.array([1], dtype=np.int64), np.array([], dtype=np.int64)]
        s1 = s2 = np.array([0, 1], dtype=np.int32)
        cuts = np.zeros((6, 1), dtype=np.float32)

        # Test partial H-Bonds and X-Bonds
        for inter in ['HBDonor', 'HBAcceptor', 'XBDonor', 'XBAcceptor']:
            others(xyz, 0, s1, s2, ball, s1, s1, s1, s1, s1, s1, s1, s1, s1,
                   s1, s1,
                   np.ones(2, dtype=np.float32), cuts, [inter], False, True,
                   np.zeros(2, dtype=np.int32))

    def test_atomic_false_branch(self):
        """Hits the same-residue filtering logic."""
        xyz = np.zeros((2, 3), dtype=np.float32)
        s1 = s2 = np.array([0, 1], dtype=np.int32)
        others(xyz, 0, s1, s2,
               [np.array([1], dtype=np.int64), np.array([], dtype=np.int64)],
               s1, s1, s1, s1, s1, s1, s1, s1, s1, s1, s1,
               np.ones(2, dtype=np.float32),
               np.zeros((6, 1), dtype=np.float32), ['VdWContact'], False,
               False, np.array([10, 10], dtype=np.int32))

    def test_others_comprehensive(self):
        np.random.seed(42)
        xyz = np.random.rand(10, 3).astype(np.float32)
        ball = [np.array([1], dtype=np.int64)] * 10
        s1 = s2 = np.arange(10, dtype=np.int32)

        # ALL interactions restored
        selected = ['VdWContact', 'CloseContact', 'Hydrophobic', 'Cationic',
                    'Anionic',
                    'MetalDonor', 'MetalAcceptor', 'HBAcceptor', 'HBDonor',
                    'XBAcceptor', 'XBDonor']
        cuts = np.full((6, len(selected)), 5.0, dtype=np.float32)

        with np.errstate(all='ignore'):
            others(xyz, 0, s1, s2, ball, s1, s1, s1, s1, s1, s1, s1, s1, s1,
                   s1, s1,
                   np.ones(10, dtype=np.float32), cuts, selected, True, True,
                   np.zeros(10, dtype=np.int32))

    def test_others_none_return(self):
        res_ijf, _ = others(None, 0, None, None, None, None, None, None, None,
                            None,
                            None, None, None, None, None, None, None, None,
                            ['None'], False, True, None)
        assert res_ijf.size == 0
