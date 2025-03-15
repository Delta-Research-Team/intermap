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
        groups = npi.group_by(ijfs[:, :2])
        sorter = groups.index.sorter
        slices = groups.index.slices
        indices = np.split(sorter, slices[1:-1])

        if not isinstance(inters, str):
            for index in indices:
                sel_ijfs = ijfs[index]
                sel_inters = inters[index]
                i, j, _ = sel_ijfs[0]
                inters_to_assign = np.any(sel_inters, axis=0).nonzero()[0]

                for inter in inters_to_assign:
                    frames = sel_ijfs[:, 2][sel_inters[:, inter]].tolist()
                    self.dict[(i, j, inter)][frames] = True
        else:
            for index in indices:
                sel_ijfs = ijfs[index]
                i, j, wat, f = sel_ijfs[0]
                frames = sel_ijfs[:, 3].tolist()
                self.dict[(i, j, wat, 'wb1')][frames] = True

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

    def save(self, path):
        """
        Save the dictionary to a csv file
        """

        def generate_lines():
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

        with open(path, 'w') as file:
            file.write('sel1_atom,sel2_atom,water_atom,interaction_name,prevalence, timeseries\n')
            for line in generate_lines():
                file.write(line)
