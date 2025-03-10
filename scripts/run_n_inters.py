# Created by rglez at 3/8/25
"""
Get No. of interactions detected in a trajectory
"""
import os
from collections import defaultdict
from os.path import join

import prolif
from bitarray.util import sc_decode
from rgpack import generals as gnl
from tqdm import tqdm


def decode_tsv(file):
    count = defaultdict(int)
    with open(file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            count[line.split('\t')[1]] += 1
    return count


def decode_pickle(file):
    imap_dict = gnl.unpickle_from_file(file)
    types = defaultdict(list)
    for key in imap_dict:
        n_frames = sum(sc_decode(imap_dict[key]['time']))
        types[key[-1]].append(n_frames)
    count = {key: sum(types[key]) for key in types}
    return count


def decode_prolif(file):
    count = defaultdict(int)
    with open(file) as f:
        for line in f:
            if line.startswith('ligand'):
                continue
            counter = 0
            for x in line.strip().split(','):
                if x == 'True':
                    counter += 1
            count[line.split(',')[2]] += counter
    return count


def get_counts(files, function):
    data = defaultdict(dict)
    for file in tqdm(files):
        print(file)
        splitted = file.split(os.sep)
        traj_name = splitted[-3]
        selection = splitted[-2]
        data[traj_name][selection] = function(file)
    return data


# =============================================================================
#
# =============================================================================
input_dir = '/media/rglez/Expansion/RoyData/intermap/benchmark/n_inters'
output_dir = '/home/rglez/RoyHub/fault/intermap/scripts'

tsv = list(gnl.recursive_finder('*tsv', input_dir))
get_contacts_data = get_counts(tsv, decode_tsv)
gnl.pickle_to_file(get_contacts_data, join(output_dir, 'gc_data.pickle'))

pickles = list(gnl.recursive_finder('*pickle', input_dir))
imap_data = get_counts(pickles, decode_pickle)
gnl.pickle_to_file(imap_data, join(output_dir, 'imap_data.pickle'))

prolifs = list(gnl.recursive_finder('interacciones.csv', input_dir))
prolif_data = get_counts(prolifs, decode_prolif)
gnl.pickle_to_file(prolif_data, join(output_dir, 'prolif_data.pickle'))
