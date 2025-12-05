import os
from collections import defaultdict

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
        df = pd.read_csv(self.imap_input, header=0, na_values=('', ' '))
        name1 = [f'{x[0]}{x[1]}' for x in df.s1.str.split('_')]
        name2 = [f'{x[0]}{x[1]}' for x in df.s2.str.split('_')]
        inter = df.inter_name.str.strip().tolist()

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


def count_diverging(diverging_dict):
    counts = defaultdict(int)
    for (a, b, c) in diverging_dict:
        counts[c] += 1
    return counts


# =============================================================================
# User-defined variables
# =============================================================================
user = os.getenv('USER')

# ----| mpro
case = 'mpro'
imap_input = f'/media/rglez/Roy5T/RoyData/intermap/1000Frames/{case}/outs_intermap/residue-12-50_InterMap.csv'
plif_input = f'/media/rglez/Roy5T/RoyData/intermap/1000Frames/{case}/prolif.pkl'
out_dir = f'/home/{user}/RoyHub/intermap/scripts/identity/{case}/'
os.makedirs(out_dir, exist_ok=True)

self = Kompare(plif_input, imap_input)
plif_dict, plif_types = self.plif2dict()
imap_dict, imap_types = self.imap2dict()

missing, diverging = self.get_missing_in_target(plif_dict, imap_dict)

gnl.pickle_to_file(plif_dict, out_dir + 'prolif_dict.pkl')
gnl.pickle_to_file(imap_dict, out_dir + 'imap_dict.pkl')

(count_diverging(plif_dict))
(count_diverging(imap_dict))
(count_diverging(diverging))

imap_keys = set(imap_dict.keys())
plif_keys = set(plif_dict.keys())

not_in_imap = sorted(set.difference(plif_keys, imap_keys))
not_in_plif = sorted(set.difference(imap_keys, plif_keys))

# =============================================================================
# Renumbering
# =============================================================================
# import prody as prd
# to_renum = '/media/rglez/Roy5T/RoyData/intermap/1000Frames/Ec_T4P/prot.pdb'
# parsed = prd.parsePDB(to_renum)
# old_resnums = parsed.getResnums()
# new_resnums = list(range(1, len(old_resnums) + 1))
# parsed.setResnums(new_resnums)
# prd.writePDB(to_renum.replace('.pdb', '_renum.pdb'), parsed)
