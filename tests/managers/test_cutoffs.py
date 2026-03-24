"""
Creator: rglez
Date: 3/22/26
Env:
Desc: Testing the CutoffsManager logic and matrix construction.
"""
from types import SimpleNamespace

import numpy as np
import pytest
from numba.typed import List

from intermap.managers.cutoffs import cutoffs, CutoffsManager, get_cutoff


# =============================================================================
# MOCKS
# =============================================================================

class MockIman:
    """
    Mocks the IndexManager (iman) behavior.
    """

    def __init__(self, inters_requested):
        self.inters_requested = inters_requested

    def get_max_vdw_dist(self):
        return 3.5


# =============================================================================
# GET_CUTOFF TESTS
# =============================================================================

class TestGetCutoff:
    def test_invalid_cutoff_name(self):
        """
        Test ValueError for bad name.
        """
        with pytest.raises(ValueError, match="is not a valid cutoff name"):
            get_cutoff("dist_cut_NonExistent")

    def test_get_from_args(self):
        """
        Triggers branch where value is pulled from user config (args).
        """
        args = {'dist_cut_Hydrophobic': 10.0}
        assert get_cutoff('dist_cut_Hydrophobic', args) == 10.0

    def test_get_from_defaults(self):
        """
        Triggers branch where value is pulled from internal constants.
        """
        assert get_cutoff('dist_cut_Hydrophobic', args={}) == cutoffs[
            'dist_cut_Hydrophobic']


# =============================================================================
# CUTOFFS MANAGER TESTS
# =============================================================================

class TestCutoffsManager:

    @pytest.fixture
    def base_args(self):
        """
        Mocking the SimpleNamespace that mimics ConfigManager.config_args
        """
        return SimpleNamespace(cutoffs={})

    def test_matrix_construction_and_filling(self, base_args):
        """
        Verifies the cuts_table (M x N matrix) is built correctly.
        Uses HBDonor to trigger dist1, dist2, minAng1, maxAng1 branches.
        Uses FaceToFace to trigger minAng2, maxAng2 branches.
        """
        iman = MockIman(['HBDonor', 'FaceToFace', 'VdWContact'])
        cm = CutoffsManager(base_args, iman)

        # cuts_internal is 6 rows (dist1, dist2, min1, max1, min2, max2)
        assert cm.cuts_internal.shape[0] == 6
        assert cm.cuts_internal.dtype == np.float32

        # Verify specific values were parsed into the matrix correctly
        hbd_idx = cm.inters_internal.index('HBDonor')
        assert cm.cuts_internal[0, hbd_idx] == cutoffs[
            'dist_cut_DA']  # distCut1
        assert cm.cuts_internal[2, hbd_idx] == cutoffs[
            'min_ang_DHA']  # minAng1

        f2f_idx = cm.inters_internal.index('FaceToFace')
        assert cm.cuts_internal[4, f2f_idx] == cutoffs[
            'min_ang_nc_FaceToFace']  # minAng2

    def test_split_logic_only_aromatic(self, base_args):
        """
        Triggers branch where non-aromatics are empty (None).
        """
        iman = MockIman(['PiStacking', 'FaceToFace'])
        cm = CutoffsManager(base_args, iman)

        # Verify Numba lists are correctly assigned 'None'
        assert list(cm.selected_others) == ['None']
        assert list(cm.selected_aro) != ['None']

        assert cm.len_others == 0
        assert cm.len_aro == 2

        assert cm.max_dist_others == 0
        assert cm.max_dist_aro > 0

    def test_split_logic_only_others(self, base_args):
        """
        Triggers branch where aromatics are empty (None).
        """
        iman = MockIman(['Hydrophobic', 'Anionic'])
        cm = CutoffsManager(base_args, iman)

        assert list(cm.selected_aro) == ['None']
        assert list(cm.selected_others) != ['None']

        assert cm.len_aro == 0
        assert cm.len_others == 2

        assert cm.max_dist_aro == 0
        assert cm.max_dist_others > 0

    def test_mixed_interactions(self, base_args):
        """
        Verifies full logic when both groups are present.
        """
        iman = MockIman(['Hydrophobic', 'PiStacking'])
        cm = CutoffsManager(base_args, iman)

        assert 'Hydrophobic' in cm.selected_others
        assert 'PiStacking' in cm.selected_aro
        assert cm.len_aro == 1
        assert cm.len_others == 1

    def test_water_bridge_coverage(self, base_args):
        """
        WaterBridge has an empty dict in parse_cutoffs.
        """
        iman = MockIman(['WaterBridge'])
        cm = CutoffsManager(base_args, iman)

        assert 'WaterBridge' in cm.selected_others
        # Verify the matrix is zeros for WaterBridge as it has no numeric cuts
        idx = cm.inters_internal.index('WaterBridge')
        assert np.all(cm.cuts_internal[:, idx] == 0)

    def test_max_dist_calculation(self, base_args):
        """
        Checks the logic for max_dist_others vs max_vdw.
        """
        iman = MockIman(['Hydrophobic'])  # Default Hydrophobic dist is 4.5
        # Mocking max_vdw as 5.0 (which is larger than 4.5)
        iman.get_max_vdw_dist = lambda: 5.0
        cm = CutoffsManager(base_args, iman)

        # max_dist_others should be max(Hydrophobic_dist, max_vdw) -> max(4.5, 5.0)
        assert cm.max_dist_others == 5.0

        # Reverse test: If max_vdw is smaller than the interaction cutoff
        iman_small = MockIman(['Hydrophobic'])
        iman_small.get_max_vdw_dist = lambda: 2.0
        cm_small = CutoffsManager(base_args, iman_small)

        # Should be max(4.5, 2.0)
        assert cm_small.max_dist_others == cutoffs['dist_cut_Hydrophobic']
