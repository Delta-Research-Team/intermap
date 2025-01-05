# Created by rglez at 12/29/24
import bitarray.util as bu


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

    def fill(self, ijfs, inters):
        """
        Fill the dictionary with the interactions

        Args:
            ijfs (ndarray): indexes of the interactions and the frames
            inters (ndarray): interactions types detected in the frames
        """
        for idx, (i, j, f) in enumerate(ijfs):
            inters_types = inters[idx].nonzero()[0]
            for inter in inters_types:
                if (i, j, inter) not in self.dict:
                    zeros = self.template.copy()
                    zeros[f] = True
                    self.dict[(i, j, inter)] = {'time': zeros}

                self.dict[(i, j, inter)]['time'][f] = True

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
            time = self.dict[key]['time']
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
