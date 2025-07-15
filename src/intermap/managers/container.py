# Created by rglez at 12/29/24

from collections import defaultdict

import numpy as np
from bitarray import util as bu
from numba import njit, types
from numba.typed import List


@njit(parallel=False, cache=True)
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


@njit(parallel=False, cache=True)
def transform_wb(ijwf):
    # Initialize containers
    to_assign = {}

    # Group the interactions
    indices = groupby_wb(ijwf)
    for index in indices.values():
        sel_ijfs = ijwf[index]
        i, j, wat, f = sel_ijfs[0]
        frames = sel_ijfs[:, 3]
        key = (i, j, wat, 'WaterBridge')
        to_assign[key] = frames
    return to_assign


@njit(parallel=False, cache=True)
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


@njit(parallel=False, cache=True)
def groupby_wb(ijfw):
    """
    Group the interactions by their indices
    """
    # Initialize typed containers
    dim1, _ = ijfw.shape
    zero32 = np.int32(0)
    pair = (zero32, zero32, zero32)
    seen = {pair}
    seen.pop()
    empty_frames = List([zero32])
    indices = {pair: empty_frames}
    del indices[pair]
    empty_frames.pop()

    for i in range(dim1):
        a, b, c, d = ijfw[i]
        pair = (a, b, c)
        if pair not in seen:
            seen.add(pair)
            indices[pair] = List([np.int32(i)])
        else:
            indices[pair].append(np.int32(i))

    # Convert the lists to numpy arrays
    key = list(indices.keys())[0]
    value = np.asarray(indices[key], dtype=np.int32)
    indices_as_array = {key: value}
    del indices_as_array[key]

    for key, value in indices.items():
        indices_as_array[key] = np.asarray(value, dtype=np.int32)

    return indices_as_array


def is_compressed(dict):
    """
    Check if the dictionary is compressed

    Args:
        dict: dict of bitarrays

    Returns:
        bool: True if the dictionary is compressed, False otherwise
    """
    seed = next(iter(dict.values()))
    try:
        bu.sc_decode(seed)
        return True
    except ValueError:
        return False


class ContainerManager:
    """
    A dictionary to store InterMap interactions
    """
    swap_inters = {
        'CloseContact': 'CloseContact',
        'VdWContact': 'VdWContact',
        'Hydrophobic': 'Hydrophobic',
        'Anionic': 'Cationic',
        'Cationic': 'Anionic',
        'MetalDonor': 'MetalAcceptor',
        'MetalAcceptor': 'MetalDonor',
        'HBAcceptor': 'HBDonor',
        'HBDonor': 'HBAcceptor',
        'XBAcceptor': 'XBDonor',
        'XBDonor': 'XBAcceptor',
        'FaceToFace': 'FaceToFace',
        'EdgeToFace': 'EdgeToFace',
        'PiStacking': 'PiStacking',
        'PiCation': 'CationPi',
        'CationPi': 'PiCation',
        'AnionPi': 'PiAnion',
        'PiAnion': 'AnionPi',
        'WaterBridge': 'WaterBridge'
    }

    def __init__(self, args, iman, cuts):

        # Parse arguments from args
        self.args = args
        self.min_prevalence = self.args.min_prevalence

        # Parse arguments from cuts
        self.cuts = cuts
        self.inter_names, self.hb_idx = self.get_inter_names()

        # Parse arguments from iman
        self.iman = iman
        self.names = (
            self.iman.resid_names if self.args.resolution == 'residue'
            else self.iman.atom_names)
        self.anotations = (
            self.iman.resid_notes if self.args.resolution == 'residue'
            else self.iman.atom_notes)
        self.waters = self.iman.waters
        self.n_frames = self.iman.traj_frames.size

        # Initialize containers
        self.detect_wb = (self.iman.waters.size > 0) and (self.hb_idx.size > 0)
        if self.detect_wb:
            self.dict = defaultdict(
                lambda: bu.sc_encode(bu.zeros(self.n_frames)))
        else:
            self.dict = defaultdict(lambda: bu.zeros(self.n_frames))

    # @profile
    def fill(self, ijfs, inters):
        """
        Fill the dictionary with the interactions

        Args:
            ijfs: np.ndarray of interacting pairs ij and their frames
            inters: np.ndarray of the interactions
        """
        if ijfs.size > 0:
            if not isinstance(inters, str):
                to_assign = transform(ijfs, inters)
            else:
                to_assign = transform_wb(ijfs)

            # Adopt a compressed representation only if detect_wb is True
            if self.detect_wb:
                for key, value in to_assign.items():
                    decoded = bu.sc_decode(self.dict[key])
                    decoded[value.tolist()] = True
                    self.dict[key] = bu.sc_encode(decoded)
            else:
                for key, value in to_assign.items():
                    self.dict[key][value.tolist()] = True

    def get_line_elements(self, dict_key):
        """
        Get the elements of an output line from the dictionary key

        Args:
            dict_key (tuple): The key of the dictionary, which can be either
                a tuple of 3 elements (s1, s2, inter) or 4 elements
                (s1, s2, wat, inter_name).

        Returns:
            s1_name (str): The name of the first selection.
            s1_note (str): The note of the first selection.
            s2_name (str): The name of the second selection.
            s2_note (str): The note of the second selection.
            s3_name (str): The name of the water selection, if applicable.
            inter_name (str): The name of the interaction.
            prevalence (float): The prevalence of the interaction.
            time (bitarray): The time series of the interaction.
        """

        # Check the length of the key and assign values accordingly
        if len(dict_key) == 3:
            s1, s2, inter = dict_key
            wat = ''
            inter_name = self.inter_names[inter]
        elif len(dict_key) == 4:
            s1, s2, wat, inter_name = dict_key
        else:
            raise ValueError('The key must have 3 or 4 elements')

        # Get the names and notes for the selections
        s1_name = self.names[s1]
        s1_note = self.anotations.get(s1, '')
        s2_name = self.names[s2]
        s2_note = self.anotations.get(s2, '')
        s3_name = self.names[wat] if wat else ''
        return s1_name, s1_note, s2_name, s2_note, s3_name, inter_name

    def rename(self):
        """
        Rename the keys of the dictionary to a readable format
        """
        shared = self.iman.shared_idx
        new_dict = {}
        compressed = is_compressed(self.dict)

        for key, value in self.dict.items():
            # Compress the time series to a bitarray if not already done
            time = bu.sc_encode(value) if not compressed else value

            # Get the selections and interaction name from the key
            s1, s2 = key[0], key[1]
            (s1_name, s1_note, s2_name, s2_note, s3_name,
             inter_name) = self.get_line_elements(key)
            standard_line = (s1_name, s1_note, s2_name, s2_note, s3_name,
                             inter_name)

            both_shared = (s1 in shared) and (s2 in shared)
            if not both_shared:
                new_dict[standard_line] = time
            else:
                swap_line = (s2_name, s2_note, s1_name, s1_note, s3_name,
                             self.swap_inters[inter_name])
                new_dict[swap_line] = time
                new_dict[standard_line] = time
        self.dict = new_dict


    def get_inter_names(self):
        """
        Get the interaction names and the indices of the hydrogen bonds

        Returns:
            inter_names: list
                The list of interaction names
            hb_idx: np.ndarray
                The indices of the hydrogen bonds
        """
        selected_others = self.cuts.selected_others
        selected_aro = self.cuts.selected_aro
        inter_names = np.asarray([x for x in selected_aro if x != 'None'] +
                                 [x for x in selected_others if x != 'None'])
        hba = np.where(inter_names == 'HBAcceptor')[0]
        hbd = np.where(inter_names == 'HBDonor')[0]
        hb_idx = np.concatenate((hba, hbd))
        return inter_names, hb_idx
