# Created by gonzalezroy at 4/14/25
from collections import defaultdict
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
from rgpack import generals as gnl

import commonplots as cmp

# =============================================================================
# Load data
# =============================================================================
user = 'gonzalezroy'
plf_path = f'/home/{user}/RoyHub/intermap/scripts/identity/prolif_dict.pkl'
imap_path = f'/home/{user}/RoyHub/intermap/scripts/identity/imap_dict.pkl'
out_dir = f'/home/{user}/RoyHub/intermap/scripts/identity/'
prolif_dict = gnl.unpickle_from_file(plf_path)
imap_dict = gnl.unpickle_from_file(imap_path)

# =============================================================================
# Gather data to plot
# =============================================================================
inters_imap = list(imap_dict.keys())
inters_prolif = list(prolif_dict.keys())
inters = sorted(set(inters_imap) | set(inters_prolif), key=lambda x: x[2])

data = defaultdict(lambda: defaultdict(list))
for inter in inters:
    if inter[-1] == 'CloseContact':
        continue
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

# %% ==========================================================================
# Plot 1: General Insights
# =============================================================================
cmp.generic_matplotlib()

fig1, ax1 = plt.subplots(dpi=300)
ax1.set_xlabel('Interaction type', fontweight='bold')
ax1.set_ylabel('Interactions detected (%)', fontweight='bold')
ax1.set_ylim(0, 105)

for i, inter in enumerate(percentage):
    both = percentage[order[i][0]]['both'][0]
    plif = percentage[order[i][0]]['plif'][0]
    imap = percentage[order[i][0]]['imap'][0]
    ax1.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
    ax1.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
    ax1.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)

ax1.grid(axis='y', lw=1, ls='--', color=cmp.c2, zorder=5)
ax1.bar(i, both, color=cmp.c3, label='Both')
ax1.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
ax1.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)

ax1.legend(loc='lower center', ncol=3, fontsize='medium', fancybox=False,
           framealpha=0.95)
ax1.set_xticks(range(len(order)))
ax1.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
plt.tight_layout()
plt.savefig(join(out_dir, 'identity-per-type-global.png'))
plt.close()

# %% ==========================================================================
# Plot 2: Detailed Insights
# =============================================================================
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax2.set_xlabel('Interaction type', fontweight='bold')
fig.text(0.0, 0.6, '# Detected interactions', va='center', rotation='vertical',
         fontweight='bold')

for i, inter in enumerate(condensed):
    both = condensed[order[i][0]]['both'][0]
    plif = condensed[order[i][0]]['plif'][0]
    imap = condensed[order[i][0]]['imap'][0]
    ax1.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
    ax2.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
    ax1.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
    ax2.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
    ax1.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)
    ax2.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)

ax1.set_ylim(0, 5500)  # outliers only
ax2.set_ylim(0, 10)  # most of the data

ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)

ax1.grid(axis='y', lw=1, ls='--', color=cmp.c2, zorder=5, alpha=cmp.alpha3)
ax2.grid(axis='y', lw=1, ls='--', color=cmp.c2, zorder=5, alpha=cmp.alpha3)

ax1.bar(i, both, color=cmp.c3, label='Both')
ax2.bar(i, both, color=cmp.c3, label='Both')
ax1.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
ax2.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
ax1.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
ax2.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
ax1.legend(loc='upper center', ncol=3, fontsize='medium', fancybox=False,
           framealpha=cmp.alpha1)

ax2.set_xticks(range(len(order)))
ax2.set_xticklabels([x[0] for x in order], rotation=45, ha='right')

d = 0.25  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color=cmp.c1, mec='k', mew=1.5, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

plt.tight_layout()
plt.savefig(join(out_dir, 'identity-per-type-detailed.png'))
plt.close()
