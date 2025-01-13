# Created by rglez at 12/29/24
import itertools as it
from collections import defaultdict
from multiprocessing import Pool
import pandas as pd
import bitarray.util as bu
import numpy as np
from scipy.sparse import csr_matrix
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

    @profile
    def fill(self, ijfs, inters):
        """
        Fill the dictionary with the interactions

        Args:
            ijfs (ndarray): indexes of the interactions and the frames
            inters (ndarray): interactions types detected in the frames
        """


        # self.dict = defaultdict(self.template.copy)
        # i, j, f = ijfs.T
        # for idx in range(ijfs.shape[0]):
        #     inters_types = inters[idx].nonzero()[0]
        #     for inter in inters_types:
        #         self.dict[(i[idx], j[idx], inter)][f[idx]] = True

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
