# Created by rglez at 1/28/25

import numpy as np

import intermap.commons as cmn


# from tests.conftest import mdtrajectory, others_containers

class TestGetIndexInter:
    def test_get_index_inter_found(self):
        """Returns correct index when interaction is found"""
        inter_name = 'Hydrophobic'
        selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
        result = cmn.indices([inter_name], selected_others)
        position = np.where(result != -1)[0][0]
        assert position == 1

    def test_get_index_inter_not_found(self):
        """Returns -1 when interaction is not found"""
        inter_name = 'Anionic'
        selected_others = np.array(['VdWContact', 'Hydrophobic', 'Cationic'])
        result = cmn.indices(inter_name, selected_others)
        position = result[result != -1]
        assert position.size == 0

    def test_get_index_inter_empty_list(self):
        """Returns -1 when selected_others list is empty"""
        inter_name = 'Hydrophobic'
        selected_others = np.array([])
        result = cmn.indices(inter_name, selected_others)
        position = result[result != -1]
        assert position.size == 0

    def test_get_index_inter_multiple_occurrences(self):
        """Returns the first index when interaction appears multiple times"""
        inter_name = 'Hydrophobic'
        selected_others = np.array(
            ['VdWContact', 'Hydrophobic', 'Hydrophobic'])
        result = cmn.indices([inter_name], selected_others)
        position = np.where(result != -1)[0][0]
        assert position == 1
