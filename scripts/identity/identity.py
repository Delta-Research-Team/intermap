from collections import defaultdict

import numpy as np
import pandas as pd
from bitarray import bitarray as ba
from rgpack import generals as gnl

# =============================================================================
# User-defined variables
# =============================================================================
prolif_pickle = '/media/gonzalezroy/Expansion/RoyData/intermap/correctness/prolif/lig-prot.pkl'
imap_csv = '/media/gonzalezroy/Expansion/RoyData/intermap/correctness/intermap/imap_lig-prot_InterMap.csv'
df = gnl.unpickle_from_file(prolif_pickle)
imap = pd.read_csv(imap_csv)

# =============================================================================
# Transforming Prolif data
# =============================================================================
prolif_dict = {}
for column in df.columns:
    r1 = column[0].split('.')[0]
    r2 = column[1].split('.')[0]
    inter = column[2]
    time = np.asarray(df[column])
    prolif_dict[(r1, r2, inter)] = time

prolif_inters = set()
for x in prolif_dict:
    r1, r2, inter = x
    prolif_inters.add(inter)

# =============================================================================
# Reducing InterMap data
# =============================================================================
imap_dict_raw = dict()
for x in imap.values:
    a1, a2, a3, inter_raw, prev, time = x
    inter = str(inter_raw)
    a1_resname, a1_resid, a1_atname = a1.split('_')
    a1_label = f'{a1_resname}{a1_resid}'
    a2_resname, a2_resid, a2_atname = a2.split('_')
    a2_label = f'{a2_resname}{a2_resid}'
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

imap_inters = set()
for x in imap_dict_raw:
    r1, r2, inter = x
    imap_inters.add(inter)

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
count = 0
for pl_inter in prolif_dict:
    if pl_inter not in imap_dict:
        count += 1
        print(f'Prolif interaction {pl_inter} not in InterMap')
print(
    f'Percentage of Prolif interactions not in InterMap: {count / len(prolif_dict) * 100:.2f}%')

# IMAP @ PLF
count = 0
for im_inter in imap_dict:
    if (im_inter not in prolif_dict) and (
            im_inter[-1] not in ['CloseContact']):
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

gnl.pickle_to_file(imap_dict, '/home/gonzalezroy/RoyHub/intermap/scripts/identity/imap_dict.pkl')
gnl.pickle_to_file(prolif_dict, '/home/gonzalezroy/RoyHub/intermap/scripts/identity/prolif_dict.pkl')
