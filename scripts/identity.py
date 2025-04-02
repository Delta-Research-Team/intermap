from collections import defaultdict

import numpy as np
import pandas as pd
from bitarray import bitarray as ba
from rgpack import generals as gnl

# import itertools as it
# pi_ring=(
#             "[a;r6]1:[a;r6]:[a;r6]:[a;r6]:[a;r6]:[a;r6]:1",
#             "[a;r5]1:[a;r5]:[a;r5]:[a;r5]:[a;r5]:1",
#         )
#
# for pi_rings in it.product(pi_ring, repeat=2):
#     print(pi_rings[0])
#     print(pi_rings[1])
#     print('---')
#             res_matches = residue.GetSubstructMatches(pi_rings[0])
#
#

# =============================================================================
# User-defined variables
# =============================================================================
prolif_pickle = '/media/rglez/Expansion/RoyData/intermap/correctness/prolif/lig-prot.pkl'
imap_pickle = '/media/rglez/Expansion/RoyData/intermap/correctness/intermap/imap_lig-prot_InterMap.csv'

df = gnl.unpickle_from_file(prolif_pickle)
imap = pd.read_csv(imap_pickle)

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
count = 0
for inter in prolif_dict:
    plf_time = prolif_dict[inter].nonzero()[0]
    imap_time = imap_dict[inter].nonzero()[0]
    equal = np.array_equal(plf_time, imap_time)
    if not equal:
        diff1 = np.setdiff1d(plf_time, imap_time)
        diff2 = np.setdiff1d(imap_time, plf_time)
        if diff1.size > 0:
            print(f'Prolif time values {diff1} not in InterMap for {inter}')
            print(plf_time)
            print(imap_time)
        if diff2.size > 0:
            count += 1
            print(f'InterMap time values {diff2} not in Prolif for {inter}')
            print(plf_time)
            print(imap_time)
# =============================================================================
# plotting
# =============================================================================
import matplotlib.pyplot as plt
import matplotlib

font = {'family': 'sans-serif',
        'size': 12}

matplotlib.rc('font', **font)
inters_imap = list(imap_dict.keys())
inters_prolif = list(prolif_dict.keys())
inters = sorted(set(inters_imap) | set(inters_prolif), key=lambda x: x[2])

data = defaultdict(lambda: defaultdict(list))
for inter in inters:
    try:
        plif_time = (prolif_dict[inter] > 0).astype(int)
    except KeyError:
        plif_time = np.zeros_like(imap_dict[inter])

    try:
        imap_time = imap_dict[inter]
    except KeyError:
        imap_time = np.zeros_like(plif_time)

    in_both = np.bitwise_and(plif_time, imap_time)
    in_plif = (plif_time - in_both).sum()
    in_imap = (imap_time - in_both).sum()

    data[inter[2]]['plif'].append(in_plif)
    data[inter[2]]['imap'].append(in_imap)
    data[inter[2]]['both'].append(in_both.sum())

condensed = defaultdict(lambda: defaultdict(list))
for inter in data:
    plif = np.array(data[inter]['plif']).sum()
    imap = np.array(data[inter]['imap']).sum()
    both = np.array(data[inter]['both']).sum()
    condensed[inter]['plif'].append(plif)
    condensed[inter]['imap'].append(imap)
    condensed[inter]['both'].append(both)

percentage = defaultdict(lambda: defaultdict(list))
for inter in condensed:
    plif = condensed[inter]['plif'][0]
    imap = condensed[inter]['imap'][0]
    both = condensed[inter]['both'][0]
    total = plif + imap + both
    percentage[inter]['plif'].append(plif / total * 100)
    percentage[inter]['imap'].append(imap / total * 100)
    percentage[inter]['both'].append(both / total * 100)

order = sorted([(x, percentage[x]['both']) for x in condensed],
               key=lambda x: x[1])

fig, ax = plt.subplots(dpi=600)
ax.set_xlabel('Interaction type', fontweight='bold')
ax.set_ylabel('Interactions detected (%)', fontweight='bold')
ax.set_ylim(0, 105)
for i, inter in enumerate(percentage):
    both = percentage[order[i][0]]['both'][0]
    plif = percentage[order[i][0]]['plif'][0]
    imap = percentage[order[i][0]]['imap'][0]
    ax.bar(i, both, color='lightgray', zorder=1, lw=0.5)
    ax.bar(i, imap, color='darkgray', bottom=both + plif, zorder=1)
    ax.bar(i, plif, color='k', bottom=both, zorder=1)

ax.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.25)
ax.bar(i, both, color='lightgray', label='Both')
ax.bar(i, imap, color='darkgray', label='InterMap only', bottom=both + plif)
ax.bar(i, plif, color='k', label='ProLIF only', bottom=both)

ax.legend(loc='upper center', ncol=3, fontsize='medium', fancybox=False,
          framealpha=0.95)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
plt.tight_layout()
plt.savefig('/home/rglez/RoyHub/intermap/scripts/identity-per-type-global.png')
plt.close()

# %%
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
# fig.subplots_adjust()

ax2.set_xlabel('Interaction type', fontweight='bold')
fig.text(0.0, 0.6, 'No. of detected interactions', va='center',
         rotation='vertical', fontweight='bold')
for i, inter in enumerate(condensed):
    both = condensed[order[i][0]]['both'][0]
    plif = condensed[order[i][0]]['plif'][0]
    imap = condensed[order[i][0]]['imap'][0]
    ax1.bar(i, both, color='lightgray', zorder=1, lw=0.5)
    ax2.bar(i, both, color='lightgray', zorder=1, lw=0.5)
    ax1.bar(i, imap, color='darkgray', bottom=both + plif, zorder=1)
    ax2.bar(i, imap, color='darkgray', bottom=both + plif, zorder=1)
    ax1.bar(i, plif, color='k', bottom=both, zorder=1)
    ax2.bar(i, plif, color='k', bottom=both, zorder=1)

ax1.set_ylim(15, 5500)  # outliers only
ax2.set_ylim(0, 15)  # most of the data

ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
# ax1.xaxis.tick_top()
ax1.tick_params(axis='both', bottom=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

ax1.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.25)
ax2.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.25)

ax1.bar(i, both, color='lightgray', label='Both')
ax2.bar(i, both, color='lightgray', label='Both')
ax1.bar(i, imap, color='darkgray', label='InterMap only', bottom=both + plif)
ax2.bar(i, imap, color='darkgray', label='InterMap only', bottom=both + plif)
ax1.bar(i, plif, color='k', label='ProLIF only', bottom=both)
ax2.bar(i, plif, color='k', label='ProLIF only', bottom=both)

ax1.legend(loc='upper center', ncol=3, fontsize='medium', fancybox=False,
           framealpha=0.85)

ax2.set_xticks(range(len(order)))
ax2.set_xticklabels([x[0] for x in order], rotation=45, ha='right')

d = 0  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1.5, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

plt.tight_layout()
plt.savefig(
    '/home/rglez/RoyHub/intermap/scripts/identity-per-type-detailed.png')
plt.close()
