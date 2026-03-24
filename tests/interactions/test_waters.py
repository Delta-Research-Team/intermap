import pytest
import numpy as np

from intermap.interactions.waters import wb1

# =============================================================================
# UNIT TESTS (Logic & Data Flow)
# =============================================================================

class TestWatersLogic:
    def test_wb1_empty_return(self):
        """Hits lines 36-37: Returns empty array if no water bridges exist."""
        # Setup data with NO waters involved
        ijf_chunk = np.array([[1, 2, 0], [3, 4, 0]], dtype=np.int32)
        inters_chunk = np.array([[True, False], [False, True]], dtype=np.bool_)
        waters = np.array([10, 11], dtype=np.int32) # Waters are indices 10, 11
        idxs = [0, 1]

        res = wb1(ijf_chunk, inters_chunk, waters, idxs, resconv=None, atomic=True)
        assert res.shape == (0, 4)

    def test_wb1_no_doubles(self):
        """Hits lines 46-55 (skipped loop): Water connects to only ONE atom."""
        # Atom 1 binds to Water 10. Water 10 binds to nothing else.
        ijf_chunk = np.array([[1, 10, 0]], dtype=np.int32)
        inters_chunk = np.array([[True, True]], dtype=np.bool_)
        waters = np.array([10], dtype=np.int32)
        idxs = [0, 1]

        res = wb1(ijf_chunk, inters_chunk, waters, idxs, resconv=None, atomic=True)
        assert res.shape == (0, 4) # No bridges formed

    def test_wb1_standard_bridge(self):
        """
        Hits lines 16-34 and 43-55: Full combination logic and column swapping.
        Oracle: Water 10 connects to Atoms 1, 2, and 3.
        Expectation: 3 possible pairs -> (1,2), (1,3), (2,3).
        """
        ijf_chunk = np.array([
            [1, 10, 0],  # Water in col 'j'
            [2, 10, 0],  # Water in col 'j'
            [10, 3, 0],  # Water in col 'i' (tests the swap block lines 29-34)
            [5, 6, 0]    # Noise (no water)
        ], dtype=np.int32)

        inters_chunk = np.array([
            [True, False],
            [False, True],
            [True, True],
            [True, False]
        ], dtype=np.bool_)

        waters = np.array([10], dtype=np.int32)
        idxs = [0, 1]

        res = wb1(ijf_chunk, inters_chunk, waters, idxs, resconv=None, atomic=True)

        # 3 atoms connected to the same water -> 3 combinations
        assert res.shape == (3, 4)

        # Column 2 should contain the water index (10), Column 3 the frame (0)
        assert np.all(res[:, 2] == 10)
        assert np.all(res[:, 3] == 0)

        # Verify the atoms connected (columns 0 and 1) are from the set {1, 2, 3}
        flat_atoms = res[:, :2].flatten()
        assert set(flat_atoms) == {1, 2, 3}

    def test_wb1_non_atomic_resolution(self):
        """Hits lines 11-13: Converts atomic indices to residue indices."""
        # When atomic=False, the input ijf_chunk contains RESIDUE indices.
        # Atom 10 belongs to Residue 100.
        ijf_chunk = np.array([
            [1, 100, 0],
            [2, 100, 0]
        ], dtype=np.int32)

        inters_chunk = np.array([[True, True], [True, False]], dtype=np.bool_)
        waters = np.array([10], dtype=np.int32) # User passed atom index 10

        resconv = np.zeros(15, dtype=np.int32)
        resconv[10] = 100 # Atom 10 -> Residue 100

        res = wb1(ijf_chunk, inters_chunk, waters, [0, 1], resconv, atomic=False)

        # Bridge formed between residues 1 and 2, mediated by residue 100
        assert res.shape == (1, 4)
        assert res[0, 2] == 100
