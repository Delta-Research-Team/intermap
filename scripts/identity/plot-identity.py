# Created by gonzalezroy at 4/14/25
import os
from collections import defaultdict
from os.path import join

import matplotlib.pyplot as plt
import numpy as np
from rgpack import generals as gnl

from scripts import commonplots as cmp


def plot_box(ax, data, positions, color):
    bplot = ax.boxplot(data, positions=positions, patch_artist=True,
                       boxprops=dict(facecolor='k', color=color,
                                     linewidth=0),
                       whiskerprops=dict(color=color),
                       capprops=dict(color=color),
                       medianprops=dict(color=color, linewidth=2),

                       flierprops=dict(marker='o', markersize=5, alpha=1,
                                       markerfacecolor=color,
                                       markeredgecolor=color),
                       showfliers=True, widths=1)
    for patch in bplot['boxes']:
        patch.set_facecolor(color)
        patch.set_alpha(0.65)

    return bplot


# =============================================================================
# Load data
# =============================================================================
user = os.getenv('USER')
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
    try:
        plif_time = set(prolif_dict[inter])
    except KeyError:
        plif_time = set()

    try:
        imap_time = set(imap_dict[inter])
    except KeyError:
        imap_time = set()

    in_both = set.intersection(plif_time, imap_time)
    in_plif = set.difference(plif_time, in_both)
    in_imap = set.difference(imap_time, in_both)

    data[inter[2]]['plif'].append(len(in_plif))
    data[inter[2]]['imap'].append(len(in_imap))
    data[inter[2]]['both'].append(len(in_both))

# %%===========================================================================
#
# =============================================================================
cmp.generic_matplotlib()

n_inters = len(data)
x_pos = range(3 * n_inters)

fig1, ax1 = plt.subplots(dpi=300, figsize=(14, 5))
ax1.set_xlabel('Interaction type', fontweight='bold')
ax1.set_ylabel('Interactions count distribution', fontweight='bold')

colors = {'both': cmp.c3, 'imap': cmp.c2, 'plif': cmp.c1}
sorted_inters = sorted(data.keys(), key=lambda x: (np.mean(data[x]['plif']),
                                                   np.mean(data[x]['both']),
                                                   np.mean(data[x]['imap'])),
                       reverse=False)
spacer = 5
current = 0
label_pos = {}
for inter in sorted_inters:
    both = data[inter]['both']
    bbplot1 = plot_box(ax1, both, [current], colors['both'])
    label_pos[inter] = current
    current += 1

    imap = data[inter]['imap']
    bbplot2 = plot_box(ax1, imap, [current], colors['imap'])
    label_pos[inter] = current
    current += 1

    plif = data[inter]['plif']
    bbplot3 = plot_box(ax1, plif, [current], colors['plif'])
    label_pos[inter] = current
    current += spacer

# Create customlegend
leg = ax1.legend(
    [bbplot1["boxes"][0], bbplot2["boxes"][0], bbplot3["boxes"][0]],
    ['Both', 'InterMap only', 'ProLIF only'],
    loc='upper center', ncol=3, fontsize='small',
    fancybox=False, framealpha=1,
    handlelength=1.5, handletextpad=0.5,
    borderpad=0.5, borderaxespad=0.5, frameon=True, bbox_to_anchor=(0.5, 1.06))

cmp.reset_matplotlib()
pos = np.asarray(list(label_pos.values())) - 1
ax1.set_xticks(pos,
               labels=list(label_pos.keys()), rotation=45, ha='right',
               fontsize=14)
ax1.set_xlim(-2, current + 1)
plt.savefig(join(out_dir, 'identity-per-type.png'),
            bbox_inches='tight')
#
# # %%===========================================================================
# #
# # =============================================================================
# condensed = defaultdict(lambda: defaultdict(list))
# for inter in data:
#     plif = len(data[inter]['plif'])
#     imap = len(data[inter]['imap'])
#     both = len(data[inter]['both'])
#     condensed[inter]['plif'].append(plif)
#     condensed[inter]['imap'].append(imap)
#     condensed[inter]['both'].append(both)
#
# percentage = defaultdict(lambda: defaultdict(list))
# for inter in condensed:
#     plif = condensed[inter]['plif'][0]
#     imap = condensed[inter]['imap'][0]
#     both = condensed[inter]['both'][0]
#     total = plif + imap + both
#     percentage[inter]['plif'].append(plif / total * 100)
#     percentage[inter]['imap'].append(imap / total * 100)
#     percentage[inter]['both'].append(both / total * 100)
#
# order = sorted([(x, percentage[x]['both']) for x in condensed],
#                key=lambda x: x[1])
#
# # %% ==========================================================================
# # Plot 1: General Insights
# # =============================================================================
# cmp.generic_matplotlib()
#
# fig1, ax1 = plt.subplots(dpi=300)
# ax1.set_xlabel('Interaction type', fontweight='bold')
# ax1.set_ylabel('Interactions detected (%)', fontweight='bold')
# ax1.set_ylim(0, 105)
#
# for i, inter in enumerate(percentage):
#     both = percentage[order[i][0]]['both'][0]
#     plif = percentage[order[i][0]]['plif'][0]
#     imap = percentage[order[i][0]]['imap'][0]
#     ax1.bar(i, both, color=cmp.c3, zorder=1, lw=0.5, alpha=cmp.alpha1)
#     ax1.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1,
#             alpha=cmp.alpha1)
#     ax1.bar(i, plif, color=cmp.c1, bottom=both, zorder=1, alpha=cmp.alpha1)
# plt.show()
# ax1.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.5)
# ax1.bar(i, both, color=cmp.c3, label='Both', alpha=0)
# ax1.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif,
#         alpha=0)
# ax1.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both, alpha=0)
#
# leg = ax1.legend(loc='lower center', ncol=3, fontsize='medium', fancybox=False,
#                  framealpha=0.95)
# for lh in leg.legend_handles:
#     lh.set_alpha(1)
# ax1.set_xticks(range(len(order)))
# ax1.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
# plt.tight_layout()
# plt.savefig(join(out_dir, 'identity-per-type-global.png'))
# plt.close()
#
# # %% ==========================================================================
# # Plot 2: Detailed Insights
# # =============================================================================
# fig, ax1 = plt.subplots(1, 1, sharex=True)
#
# # ax2.set_xlabel('Interaction type', fontweight='bold')
# fig.text(0.0, 0.6, '# Detected interactions', va='center', rotation='vertical',
#          fontweight='bold')
#
# for i, inter in enumerate(condensed):
#     both = condensed[order[i][0]]['both'][0]
#     plif = condensed[order[i][0]]['plif'][0]
#     imap = condensed[order[i][0]]['imap'][0]
#     ax1.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
#     # ax2.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
#     ax1.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
#     # ax2.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
#     ax1.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)
#     # ax2.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)
#
# # ax1.set_ylim(0, 100)  # outliers only
# # ax2.set_ylim(0, 10)  # most of the data
#
# ax1.spines.bottom.set_visible(False)
# # ax2.spines.top.set_visible(False)
#
# ax1.grid(axis='y', lw=1, ls='--', color='gray', zorder=5, alpha=cmp.alpha3)
# # ax2.grid(axis='y', lw=1, ls='--', color=cmp.c2, zorder=5, alpha=cmp.alpha3)
#
# ax1.bar(i, both, color=cmp.c3, label='Both')
# # ax2.bar(i, both, color=cmp.c3, label='Both')
# ax1.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
# # ax2.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
# ax1.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
# # ax2.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
# ax1.legend(loc='upper center', ncol=3, fontsize='medium', fancybox=False,
#            framealpha=cmp.alpha1)
#
# # ax2.set_xticks(range(len(order)))
# # ax2.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
#
# d = 0.25  # proportion of vertical to horizontal extent of the slanted line
# kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
#               linestyle="none", color=cmp.c1, mec='k', mew=1.5, clip_on=False)
# ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
# # ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
#
# plt.tight_layout()
# plt.savefig(join(out_dir, 'identity-per-type-detailed.png'))
# plt.close()
