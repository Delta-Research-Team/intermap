import numpy as np
import pytest
from bitarray import util as bu

from intermap.managers.container import (ContainerManager, groupby, groupby_wb,
                                         is_compressed, transform,
                                         transform_wb)


# =============================================================================
# NUMBA KERNEL TESTS (Data Grouping & Transformation)
# =============================================================================

class TestContainerTransforms:
    def test_groupby_and_transform(self):
        """Tests the core Numba grouping logic."""
        # Note: Numba explicitly requires int64 for the ijfs array signature
        ijfs = np.array([[1, 2, 0], [1, 2, 1], [3, 4, 0]], dtype=np.int64)
        inters = np.array([[True, False], [False, True], [True, True]],
                          dtype=np.bool_)

        # Test groupby
        idx_dict = groupby(ijfs)
        assert len(idx_dict) == 2
        assert np.array_equal(idx_dict[(1, 2)], [0, 1])

        # Test transform
        res = transform(ijfs, inters)
        assert list(res[(1, 2, 0)]) == [0]  # Pair 1-2, Inter 0 -> Frame 0
        assert list(res[(1, 2, 1)]) == [1]  # Pair 1-2, Inter 1 -> Frame 1
        assert list(res[(3, 4, 0)]) == [0]  # Pair 3-4, Inter 0 -> Frame 0

    def test_groupby_wb_and_transform_wb(self):
        """Tests the WaterBridge Numba grouping logic."""
        ijfw = np.array([[1, 2, 99, 0], [1, 2, 99, 1]], dtype=np.int64)

        # Test groupby_wb
        idx_dict = groupby_wb(ijfw)
        assert len(idx_dict) == 1
        assert np.array_equal(idx_dict[(1, 2, 99)], [0, 1])

        # Test transform_wb
        res = transform_wb(ijfw)
        assert list(res[(1, 2, 99, 'WaterBridge')]) == [0, 1]


class TestContainerCompression:
    def test_is_compressed(self):
        """Hits the try/except block for bitarray compression."""
        # Compressed: Create a valid bitarray and encode it properly
        valid_b = bu.zeros(10)
        comp_dict = {'k': bu.sc_encode(valid_b)}
        assert is_compressed(comp_dict) == True

        # Uncompressed (Forces ValueError by using invalid sc_encode format)
        # sc_decode strictly expects a specific byte format.
        # Passing an ordinary, un-encoded bitarray of all 1s forces it to crash.
        invalid_b = bu.bitarray('1111111111')
        uncomp_dict = {'k': invalid_b}
        assert is_compressed(uncomp_dict) == False


# =============================================================================
# CONTAINER MANAGER TESTS
# =============================================================================

class MockArgs:
    min_prevalence = 0.0
    resolution = 'atom'


class MockCuts:
    selected_aro = ['FaceToFace', 'None']
    selected_others = ['HBDonor', 'None']


class MockIman:
    """Provides surgical control over naming dictionaries."""
    resid_names = {1: 'RES1', 2: 'RES2'}
    atom_names = {1: 'ATOM1', 2: 'ATOM2', 99: 'WAT'}
    resid_notes = {1: 'note1', 2: ''}
    atom_notes = {1: 'anote1', 2: ''}
    waters = np.array([99])
    traj_frames = np.arange(3)
    shared_idx = np.array([1, 2])  # Forces Both_Shared=True in rename()


class TestContainerManager:
    def test_full_flow_and_rename_swap(self):
        """
        Tests initialization, filling, and the complex rename logic.
        Hits both_shared=True branch.
        """
        cm = ContainerManager(MockArgs(), MockIman(), MockCuts())

        # inter_names = ['FaceToFace', 'HBDonor']
        # Indices: FaceToFace=0, HBDonor=1

        # 1. Fill standard interactions
        ijfs = np.array([[1, 2, 0], [1, 2, 1]], dtype=np.int64)
        inters = np.array([[True, False], [False, True]], dtype=np.bool_)
        cm.fill(ijfs, inters)

        # 2. Fill water bridges
        ijfw = np.array([[1, 2, 99, 2]], dtype=np.int64)
        cm.fill(ijfw, 'wb')

        # Verify un-renamed keys exist
        assert (1, 2, 0) in cm.dict
        assert (1, 2, 99, 'WaterBridge') in cm.dict

        # 3. Test line element parser
        names3 = cm.get_line_elements((1, 2, 0))
        assert names3 == ('ATOM1', 'anote1', 'ATOM2', '', '', 'FaceToFace')

        names4 = cm.get_line_elements((1, 2, 99, 'WaterBridge'))
        assert names4 == ('ATOM1', 'anote1', 'ATOM2', '', 'WAT', 'WaterBridge')

        with pytest.raises(ValueError, match="key must have 3 or 4 elements"):
            cm.get_line_elements((1, 2, 3, 4, 5))

        # 4. Test Rename (Swap Logic)
        cm.rename()

        # Since 1 and 2 are in shared_idx, HBDonor creates a duplicate swapped line for HBAcceptor
        expected_std_hbd = ('ATOM1', 'anote1', 'ATOM2', '', '', 'HBDonor')
        expected_swp_hba = ('ATOM2', '', 'ATOM1', 'anote1', '', 'HBAcceptor')

        assert expected_std_hbd in cm.dict
        assert expected_swp_hba in cm.dict

    def test_uncompressed_and_empty_fill(self):
        """Hits detect_wb=False and empty array branches."""
        iman_no_wat = MockIman()
        iman_no_wat.waters = np.array([])

        cm = ContainerManager(MockArgs(), iman_no_wat, MockCuts())
        assert cm.detect_wb == False

        # Hit empty fill
        cm.fill(np.array([], dtype=np.int64), np.array([]))
        assert len(cm.dict) == 0

        # Hit standard fill uncompressed
        cm.fill(np.array([[1, 2, 0]], dtype=np.int64),
                np.array([[True, False]], dtype=np.bool_))
        cm.rename()
        assert len(cm.dict) > 0

    def test_residue_resolution(self):
        """Hits resolution='residue' attribute assignment."""
        args = MockArgs()
        args.resolution = 'residue'
        cm = ContainerManager(args, MockIman(), MockCuts())

        assert cm.names == MockIman.resid_names
        assert cm.anotations == MockIman.resid_notes

    def test_uncompressed_fill_and_get_names(self):
        """Hits uncompressed 'else' block and 'get_inter_names'."""
        # 1. Setup Mock where detect_wb is False
        iman_no_wat = MockIman()
        iman_no_wat.waters = np.array([])

        cuts = MockCuts()
        cuts.selected_aro = ['FaceToFace']
        cuts.selected_others = ['CloseContact']

        cm = ContainerManager(MockArgs(), iman_no_wat, cuts)

        # 2. Fill with uncompressed data (Hits line 141-142)
        ijfs = np.array([[1, 2, 0]], dtype=np.int64)
        inters = np.array([[True]], dtype=np.bool_)
        cm.fill(ijfs, inters)

        # Verify it was added
        assert (1, 2, 0) in cm.dict

        # 3. Test get_inter_names logic (Hits line 284)
        names, hb_idx = cm.get_inter_names()
        assert 'FaceToFace' in names
        assert 'CloseContact' in names
        assert len(hb_idx) == 0  # No H-bonds selected
