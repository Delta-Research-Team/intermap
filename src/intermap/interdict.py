# Created by rglez at 12/29/24
import bitarray.util as bu


class InterDict:
    """
    A dictionary to store InterMap interactions
    """

    def __init__(self, num_frames):
        self.n_frames = num_frames
        self.template = bu.zeros(num_frames)
        self.dict = {}

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

    def get_prevalence(self):
        """
        Get the prevalence of the interactions
        """
        prevalence = {}
        for key in self.dict:
            count = self.dict[key]['time'].count()
            percentage = round(count / self.n_frames * 100, 2)
            self.dict[key]['prevalence'] = percentage
        return prevalence

    def compress(self):
        """
        Compress the bit dictionary
        """
        for key in self.dict:
            self.dict[key]['time'] = bu.sc_encode(self.dict[key]['time'])
