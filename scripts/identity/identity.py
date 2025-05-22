import numpy as np
import pandas as pd
from bitarray import bitarray as ba
from rgpack import generals as gnl


class Kompare:
    def __init__(self, plif_pickle, imap_csv, ignores=()):
        self.plif_input = plif_pickle
        self.imap_input = imap_csv
        self.ignores = ignores

    def plif2dict(self):
        """
        Convert a Prolif pickle file to a dictionary.
        """
        df = gnl.unpickle_from_file(self.plif_input)
        prolif_dict = {}
        for column in df.columns:
            r1 = column[0].split('.')[0].strip()
            r2 = column[1].split('.')[0].strip()
            inter = column[2].strip()
            time = np.asarray(df[column], dtype=bool).nonzero()[0]
            if r1 == r2:
                continue
            prolif_dict[(r1, r2, inter)] = time
        inters_types = set([x[-1] for x in prolif_dict.keys()])
        return prolif_dict, inters_types

    def imap2dict(self):
        """
        Convert an InterMap csv file to a dictionary.
        """
        df = pd.read_csv(self.imap_input, header=2, na_values=('', ' '))
        name1 = [f'{x[0]}{x[1]}' for x in df.sel1.str.split('_')]
        name2 = [f'{x[0]}{x[1]}' for x in df.sel2.str.split('_')]
        inter = df.interaction_name.str.strip().tolist()

        imap_dict = {}
        for i, x in enumerate(name1):
            key = (name1[i].strip(), name2[i].strip(), inter[i].strip())
            time = np.asarray(list(ba(df.iloc[i, -1]))).nonzero()[0]
            imap_dict[key] = time
        inters_types = set([x[-1] for x in imap_dict.keys()])
        return imap_dict, inters_types

    def get_missing_in_target(self, reference, target):
        """
        Get the missing interactions in the target dictionary compared to the reference.
        """
        missing = {}
        diverging = {}
        for inter in reference:
            type_ = inter[-1]
            if type_ in self.ignores:
                continue
            ref_time = reference.get(inter)
            tar_time = target.get(inter, None)

            if tar_time is None:
                missing[inter] = ref_time
            else:
                equal = np.array_equal(ref_time, tar_time)
                if not equal:
                    difference = np.setdiff1d(ref_time, tar_time)
                    diverging[inter] = difference
        return missing, diverging


# =============================================================================
# User-defined variables
# =============================================================================
out_dir = '/home/gonzalezroy/RoyHub/intermap/scripts/identity/'
imap_input = '/media/gonzalezroy/Roy2TB/RoyData/intermap/IDENTITY-FINAL/imaps/mpro/mpro_InterMap_full.csv'
plif_input = '/media/gonzalezroy/Roy2TB/RoyData/intermap/IDENTITY-FINAL/mpro/prolif.pkl'

kompare = Kompare(plif_input, imap_input)
plif_dict, plif_types = kompare.plif2dict()
imap_dict, imap_types = kompare.imap2dict()
missing, diverging = kompare.get_missing_in_target(plif_dict, imap_dict)

gnl.pickle_to_file(plif_dict, out_dir + 'prolif_dict.pkl')
gnl.pickle_to_file(imap_dict, out_dir + 'imap_dict.pkl')
