# Created by rglez at 12/29/24
from collections import defaultdict

import bitarray.util as bu
import numpy as np
import numpy_indexed as npi


# todo: numba the fill method of InterDict
class InterDict:
    """
    A dictionary to store InterMap interactions
    """

    def __init__(self, format_, min_prevalence, frames_id, atom_names,
                 inter_names):

        # Parse arguments
        self.format = format_
        self.min_prevalence = min_prevalence
        self.atom_names = atom_names
        self.inter_names = inter_names
        self.n_frames = frames_id.size

        # Initialize containers
        self.dict = {}
        self.template = bu.zeros(self.n_frames)

    # @profile
    def fill(self, ijfs, inters):
        groups = npi.group_by(ijfs[:, :2])
        sorter = groups.index.sorter
        slices = groups.index.slices

        self.dict = defaultdict(self.template.copy)
        indices = iter(np.split(sorter, slices[1:-1]))
        for index in indices:
            sel_ijfs = ijfs[index]
            sel_inters = inters[index]
            i, j, _ = sel_ijfs[0]

            inters_to_assign = np.any(sel_inters, axis=0).nonzero()[0]
            for inter in inters_to_assign:
                frames = sel_ijfs[:, 2][sel_inters[:, inter]].tolist()
                self.dict[(i, j, inter)][frames] = True

    def pack(self):
        """
        Get the prevalence of the interactions
        """

        packed_dict = {}
        keys = list(self.dict.keys())
        for key in keys:
            s1_id, s2_id, inter_id = key
            s1_name = self.atom_names[s1_id]
            s2_name = self.atom_names[s2_id]
            inter_name = self.inter_names[inter_id]
            time = self.dict[key]
            prevalence = round(time.count() / self.n_frames * 100, 2)
            del self.dict[key]

            if prevalence < self.min_prevalence:
                continue

            if self.format == 'extended':
                packed = bu.sc_encode(time)
                packed_dict[(s1_name, s2_name, inter_name)] = {
                    'time': packed, 'prevalence': prevalence}
            else:
                packed_dict[(s1_name, s2_name, inter_name)] = prevalence

        self.dict = packed_dict
