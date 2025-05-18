import os
from os.path import basename, join

import numpy as np
import pandas as pd
from bitarray import bitarray as ba
from rgpack import generals as gnl

# =============================================================================
# User-defined variables
# =============================================================================
root_dir = '/media/rglez/Roy2TB/Dropbox/DockingBoys/storage/identity_prolif_100frames/'
cases = [join(root_dir, x) for x in os.listdir(root_dir)]

results = {}
for case in cases:
    pkl = next(gnl.recursive_finder('fingerprint.pkl', case))
    csv = next(gnl.recursive_finder('*_full.csv', case))
    results[basename(case)] = (pkl, csv)

case = 'mpro'

df = gnl.unpickle_from_file(results[case][0]).to_dataframe()
imap = pd.read_csv(results[case][1], header=1, na_values=('', ' '))

# =============================================================================
# Transforming Prolif data
# =============================================================================
prolif_dict = {}
for column in df.columns:
    r1 = column[0].split('.')[0]
    r2 = column[1].split('.')[0]
    inter = column[2]
    time = np.asarray(df[column], dtype=bool)
    prolif_dict[(r1, r2, inter)] = time
total_prolif = 0
for x in prolif_dict:
    total_prolif += prolif_dict[x].sum()

prolif_inters = set()
for x in prolif_dict:
    r1, r2, inter = x
    prolif_inters.add(inter.strip())

# =============================================================================
# Reducing InterMap data
# =============================================================================
imap_dict_raw = dict()
for x in imap.values:
    a1, _, a2, _, a3, inter_raw, prev, time = x
    a1, a2, inter_raw = a1.strip(), a2.strip(), inter_raw.strip()
    inter = str(inter_raw).strip()
    a1_resname, a1_resnum, a1_resid = a1.split('_')
    a1_label = f'{a1_resname}{a1_resnum}'
    a2_resname, a2_resnum, a2_resid = a2.split('_')
    a2_label = f'{a2_resname}{a2_resnum}'.strip()
    current_time = ba(time)

    if (a1_label, a2_label, inter) in imap_dict_raw:
        previous_time = imap_dict_raw[(a1_label, a2_label, inter)]
        or_time = previous_time | current_time
        imap_dict_raw[(a1_label, a2_label, inter)] = or_time
    else:
        imap_dict_raw[(a1_label, a2_label, inter)] = current_time

imap_dict = dict()
for x in imap_dict_raw:
    imap_dict[x] = np.asarray(imap_dict_raw[x].tolist())
total_imap = 0
for x in imap_dict:
    total_imap += imap_dict[x].sum()

imap_inters = set()
for x in imap_dict_raw:
    r1, r2, inter = x
    imap_inters.add(inter.strip())

# =============================================================================
# Comparing Prolif and InterMap data in general terms
# =============================================================================

# Interactions types in Prolif but not in InterMap
not_imap = prolif_inters - imap_inters
print('Interactions in Prolif but not in InterMap:', not_imap)

# Interactions in InterMap but not in Prolif
not_prolif = imap_inters - prolif_inters
print('Interactions in InterMap but not in Prolif:', not_prolif)

# PLF @ IMAP
# todo: what about self-interactions? within same residue? Remove a==b !
INTER = 'CationPi'
count = 0
for pl_inter in prolif_dict:
    a, b, c = pl_inter
    pl_inverse = (b, a, c)
    if (c != INTER) or (a == b):
        continue
    if pl_inter not in imap_dict:
        count += 1
        print(f'Prolif interaction {pl_inter} not in InterMap')
print(
    f'Percentage of Prolif {INTER} not in InterMap: {count / len(prolif_dict) * 100:.2f}%')

# IMAP @ PLF
count = 0
for im_inter in imap_dict:
    d, e, f = im_inter
    im_inverse = (e, d, f)
    if (d == e) or (f != INTER):
        continue
    if (im_inter not in prolif_dict) and (
            im_inter[-1] in ['CloseContact']):
        count += 1
        print(f'InterMap interaction {im_inter} not in Prolif')
print(
    f'Percentage of InterMap interactions not in Prolif: {count / len(imap_dict) * 100:.2f}%')

# =============================================================================
# Comparing Prolif and InterMap data in specific terms
# =============================================================================
# count = 0
# for inter in prolif_dict:
#     plf_time = prolif_dict[inter].nonzero()[0]
#     imap_time = imap_dict[inter].nonzero()[0]
#     equal = np.array_equal(plf_time, imap_time)
#     if not equal:
#         diff1 = np.setdiff1d(plf_time, imap_time)
#         diff2 = np.setdiff1d(imap_time, plf_time)
#         if diff1.size > 0:
#             print(f'Prolif time values {diff1} not in InterMap for {inter}')
#             print(plf_time)
#             print(imap_time)
#         if diff2.size > 0:
#             count += 1
#             print(f'InterMap time values {diff2} not in Prolif for {inter}')
#             print(plf_time)
#             print(imap_time)

gnl.pickle_to_file(imap_dict,
                   '/home/gonzalezroy/RoyHub/intermap/scripts/identity/imap_dict.pkl')
gnl.pickle_to_file(prolif_dict,
                   '/home/gonzalezroy/RoyHub/intermap/scripts/identity/prolif_dict.pkl')
