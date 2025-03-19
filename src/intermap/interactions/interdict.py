# Created by rglez at 12/29/24

from collections import defaultdict

import bitarray.util as bu
import numpy as np
from numba import njit, types
from numba.typed import List


@njit(parallel=False, cache=False)
def transform(ijfs, inters):
    # Initialize containers
    i = 0
    j = 0
    inter = 0
    key = (i, j, inter)
    value = np.asarray([i, j, inter], dtype=np.int32)
    to_assign = {key: value}
    del to_assign[key]

    # Group the interactions
    indices = groupby(ijfs)
    for index in indices.values():
        sel_ijfs = ijfs[index]
        sel_inters = inters[index]
        i, j, _ = sel_ijfs[0]
        inters_to_assign = sel_inters.sum(axis=0).nonzero()[0]
        for inter in inters_to_assign:
            frames = sel_ijfs[:, 2][sel_inters[:, inter]]
            key = (i, j, inter)
            to_assign.update({key: frames})
    return to_assign


@njit(parallel=False, cache=False)
def transform_wb(ijfs):
    # Initialize containers
    i = 0
    j = 0
    inter = 0
    key = (i, j, inter)
    value = np.asarray([i, j, inter], dtype=np.int32)
    to_assign = {key: value}
    del to_assign[key]

    # Group the interactions
    indices = groupby(ijfs)
    for index in indices.values():
        sel_ijfs = ijfs[index]
        i, j, wat, f = sel_ijfs[0]
        frames = sel_ijfs[:, 3]
        key = (i, j, wat, 'wb1')
        to_assign.update({key: frames})
    return to_assign


@njit(parallel=False, cache=False)
def groupby(ijfs: types.Array(types.int64, 2, 'C')) -> types.DictType(
    types.UniTuple(types.int64, 2), types.ListType(types.int32)):
    """
    Group the interactions by their indices
    """
    # Initialize containers
    a, b = 0, 0
    pair = (a, b)
    seen = {pair}
    seen.pop()
    empty_frames = List([a])
    empty_frames.pop()
    indices = {pair: empty_frames}
    del indices[pair]

    # Group the interactions
    for i in range(ijfs.shape[0]):
        a, b, c = ijfs[i]
        pair = (a, b)
        if pair not in seen:
            seen.add(pair)
            indices[pair] = List([i])
        else:
            indices[pair].append(i)

    # Convert the lists to numpy arrays
    key = list(indices.keys())[0]
    value = np.asarray(indices[key], dtype=np.int32)
    indices_as_array = {key: value}
    del indices_as_array[key]

    for key, value in indices.items():
        indices_as_array[key] = np.asarray(value, dtype=np.int32)

    return indices_as_array


# todo: numba the fill method of InterDict
class InterDict:
    """
    A dictionary to store InterMap interactions
    """

    def __init__(self, format_, min_prevalence, frames_id, atom_names,
                 inter_names, resids, waters):

        # Parse arguments
        self.format = format_
        self.min_prevalence = min_prevalence
        self.atom_names = atom_names
        self.inter_names = inter_names
        self.n_frames = frames_id.size
        self.resids = resids
        self.waters = waters

        # Initialize containers
        self.dict = defaultdict(lambda: bu.zeros(self.n_frames))

    def fill(self, ijfs, inters):
        if not isinstance(inters, str):
            to_assign = transform(ijfs, inters)
        else:
            to_assign = transform_wb(ijfs)
        for key, value in to_assign.items():
            self.dict[key][value.tolist()] = True

    def lock(self):
        """
        Lock the dictionary to avoid further modifications
        """
        self.dict = dict(self.dict)

    def pack(self):
        """
        Get the prevalence of the interactions
        """

        packed_dict = {}
        keys = list(self.dict.keys())
        for key in keys:
            if len(key) == 3:
                s1_id, s2_id, inter_id = key
                s1_name = self.atom_names[s1_id]
                s2_name = self.atom_names[s2_id]
                inter_name = self.inter_names[inter_id]
                key2save = (s1_name, s2_name, inter_name)

            elif len(key) == 4:
                s1_id, s2_id, wat_id, inter_id = key
                s1_name = self.atom_names[s1_id]
                s2_name = self.atom_names[s2_id]
                wat_name = self.atom_names[wat_id]
                key2save = (s1_name, s2_name, wat_name, inter_id)
            else:
                raise ValueError('The key must have 3 or 4 elements')

            time = self.dict[key]
            prevalence = round(time.count() / self.n_frames * 100, 2)
            del self.dict[key]
            if prevalence < self.min_prevalence:
                continue

            if self.format == 'extended':
                packed_dict[key2save] = {'time': bu.sc_encode(time),
                                         'prevalence': prevalence}
            else:
                packed_dict[key2save] = prevalence

        self.dict = packed_dict

    def generate_lines(self):
        for key in self.dict:
            if len(key) == 3:
                s1, s2, inter = key
                wat = ''
                inter_name = self.inter_names[inter]
            elif len(key) == 4:
                s1, s2, wat, inter_name = key
            else:
                raise ValueError('The key must have 3 or 4 elements')

                s1_name = self.atom_names[s1]
                s2_name = self.atom_names[s2]
                s3_name = self.atom_names[wat] if wat else ''
                time = self.dict[key]
                prevalence = round(time.count() / self.n_frames * 100, 2)

            # Yield each line as a generator
            yield f'{s1_name},{s2_name},{s3_name},{inter_name},{prevalence},{time.to01()}\n'

    def save(self, path):
        """
        Save the dictionary to a csv file
        """

        with open(path, 'w') as file:
            file.write(
                'sel1_atom,sel2_atom,water_atom,interaction_name,prevalence, timeseries\n')
            for line in self.generate_lines():
                file.write(line)
