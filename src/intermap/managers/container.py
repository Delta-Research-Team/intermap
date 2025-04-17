# Created by rglez at 12/29/24

from collections import defaultdict
from os.path import basename

import bitarray.util as bu
import numpy as np
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


class ContainerManager:
    """
    A dictionary to store InterMap interactions
    """

    def __init__(self, args, iman, cuts):

        # Parse arguments from args
        self.args = args
        self.format = self.args.format
        self.min_prevalence = self.args.min_prevalence

        # Parse arguments from cuts
        self.cuts = cuts
        self.inter_names, self.hb_idx = self.get_inter_names()

        # Parse arguments from iman
        self.iman = iman
        self.names = (
            self.iman.resid_names if self.args.resolution == 'residue'
            else self.iman.atom_names)
        self.waters = self.iman.waters
        self.n_frames = self.iman.traj_frames.size

        # Initialize containers
        self.dict = defaultdict(lambda: bu.zeros(self.n_frames))

        # Detect water bridges ?
        self.detect_wb = (self.iman.waters.size > 0) and (self.hb_idx.size > 0)

    # @profile
    def fill(self, ijfs, inters):
        if ijfs.size > 0:
            if not isinstance(inters, str):
                to_assign = transform(ijfs, inters)
            else:
                to_assign = transform_wb(ijfs)
            for key, value in to_assign.items():
                self.dict[key][value.tolist()] = True

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

            s1_name = self.names[s1]
            s2_name = self.names[s2]
            s3_name = self.names[wat] if wat else ''
            time = self.dict[key]
            prevalence = round(time.count() / self.n_frames * 100, 2)

            # Yield each line as a generator
            yield f'{s1_name},{s2_name},{s3_name},{inter_name},{prevalence},{time.to01()}\n'

    def save(self, path):
        """
        Save the dictionary to a csv file
        """

        with open(path, 'w') as file:
            file.write(f'# {basename(self.args.topology)}\n')
            file.write(
                'sel1_atom,sel2_atom,water_atom,interaction_name,prevalence, timeseries\n')
            for line in self.generate_lines():
                file.write(line)

    def get_inter_names(self):
        selected_others = self.cuts.selected_others
        selected_aro = self.cuts.selected_aro
        inter_names = np.asarray([x for x in selected_aro if x != 'None'] +
                                 [x for x in selected_others if x != 'None'])
        hba = np.where(inter_names == 'HBAcceptor')[0]
        hbd = np.where(inter_names == 'HBDonor')[0]
        hb_idx = np.concatenate((hba, hbd))
        return inter_names, hb_idx

    # def to_graph(self):
    #     """
    #     Convert the dictionary to a graph
    #     """
    #     # Create a dictionary with the average prevalence
    #     shrinked = defaultdict(lambda: bu.zeros(self.n_frames))
    #     for key in self.dict:
    #         r1, r2 = key[0], key[1]
    #         time = self.dict[key]
    #         shrinked[(r1, r2)] |= time
    #
    #     shrinked = {key: round(value.count() / len(value), 2) * 100
    #                 for key, value in shrinked.items()}
    #
    #     # Count nodes
    #     nodes = []
    #     for r1, r2 in shrinked.keys():
    #         nodes.append(r1)
    #         nodes.append(r2)
    #     nodes_count = Counter(nodes)
    #
    #     # Create a graph from the dictionary
    #     G = nx.Graph()
    #     for key, value in shrinked.items():
    #         r1, r2 = key
    #         if value >= self.min_prevalence:
    #             # Add nodes with their degree
    #             if not G.has_node(r1):
    #                 G.add_node(r1, name=self.names[r1], degree=nodes_count[r1])
    #             if not G.has_node(r2):
    #                 G.add_node(r2, name=self.names[r2], degree=nodes_count[r2])
    #             # Add edges with their weight
    #             G.add_edge(r1, r2, weight=value)
    #
    #     # Convert to cytoscape
    #     s1_resids = self.iman.universe.select_atoms(self.iman.sel1).resindices
    #     out_name = join(self.args.output_dir, "InterGraph.csv")
    #     with open(out_name, 'w') as file:
    #         file.write('source, target, r1_sel, r2_sel, source_count, target_count, weight\n')
    #         for u, v, data in G.edges(data=True):
    #             source_count = nodes_count[u]
    #             target_count = nodes_count[v]
    #             r1_sel = 's1' if u in s1_resids else 's2'
    #             r2_sel = 's1' if v in s1_resids else 's2'
    #             file.write(f'{u},{v},{r1_sel},{r2_sel},{source_count},{target_count},{data["weight"]}\n')
    #
    #     return G

    def pack(self):
        """
        Destructive packing of the dictionary
        """
        self.dict = {key: bu.sc_encode(self.dict[key]) for key in self.dict}
