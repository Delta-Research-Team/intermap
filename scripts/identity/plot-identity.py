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

                       flierprops=dict(marker='_', markersize=5, alpha=1,
                                       markerfacecolor=color,
                                       markeredgecolor=color),
                       showfliers=True, widths=1,
                       zorder=5)
    for patch in bplot['boxes']:
        patch.set_facecolor(color)
        patch.set_alpha(1)

    return bplot


def scatter_points(ax, data, x_center, color, alpha=0.1, size=15,
                   spread=0.75):
    """
    Draw jittered scatter points around a boxplot location.

    Parameters
    ----------
    ax : matplotlib axis
    data : list of numbers
    x_center : float
        X position of the boxplot
    color : matplotlib color
    alpha : float
        Transparency
    size : float
        Marker size
    spread : float
        Maximum horizontal jitter (+/-)
    """
    xs = [x_center + np.random.uniform(-spread, spread) for _ in data]
    ax.scatter(xs, data, color=color, alpha=alpha, s=size, edgecolors='none',
               zorder=2)


# =============================================================================
# Load data
# =============================================================================
user = os.getenv('USER')

case = 'mpro'
plf_path = f'/home/{user}/RoyHub/intermap/scripts/identity/{case}/prolif_dict.pkl'
imap_path = f'/home/rglez/RoyHub/intermap/scripts/identity/{case}/imap_dict.pkl'
out_dir = f'/home/{user}/RoyHub/intermap/scripts/identity/{case}'
prolif_dict = gnl.unpickle_from_file(plf_path)
imap_dict = gnl.unpickle_from_file(imap_path)

# =============================================================================
# Gather data to plot
# =============================================================================
inters_imap = list(imap_dict.keys())
inters_prolif = list(prolif_dict.keys())
inters = sorted(set(inters_imap) | set(inters_prolif), key=lambda x: x[2])

inters_fig = set([x[-1] for x in inters_imap])

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


# ===========================================================================
# Two top panels with broken y-axis for counts (LOW_MAX = 250)
# ===========================================================================

cmp.generic_matplotlib()

import matplotlib as mpl
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['mathtext.fontset'] = 'stix'

ALPHA = 1   # uniform transparency everywhere
GRID_ALPHA = 0.15
SCATTER_ALPHA = 0.35
BOX_ALPHA = 0.55
BAR_ALPHA = 0.55

LOW_MAX = 250  # max value for lower top panel

colors = {'both': cmp.c3, 'imap': cmp.c2, 'plif': cmp.c1}
sorted_inters = [
    'CationPi', 'PiCation', 'Anionic', 'Cationic',
    'VdWContact', 'Hydrophobic',
    'HBAcceptor', 'HBDonor',
    'EdgeToFace', 'FaceToFace'
]

fig1 = plt.figure(figsize=(9, 9), dpi=300)
gs = fig1.add_gridspec(
    nrows=3,
    height_ratios=[0.25, 0.25, 0.50],
    hspace=0.1
)

ax_top_high = fig1.add_subplot(gs[0])
ax_top_low = fig1.add_subplot(gs[1], sharex=ax_top_high)
ax_main = fig1.add_subplot(gs[2], sharex=ax_top_high)

# Clean x-axis for top panels
ax_top_high.tick_params(axis='x', labelbottom=False)
ax_top_low.tick_params(axis='x', labelbottom=False)

# Labels for main panel
ax_main.set_xlabel('Interaction type')
ax_main.set_ylabel('# Frames')

# Light gridlines and spines
for ax in [ax_top_high, ax_top_low, ax_main]:
    ax.grid(axis='y', linestyle='--', linewidth=1,
            color='k', alpha=GRID_ALPHA)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# Broken axis look: hide touching spines
ax_top_high.spines['bottom'].set_visible(False)
ax_top_low.spines['top'].set_visible(False)

spacer = 5
current = 0
label_pos = {}

# Collect bar data first
bar_positions = []
bar_counts = []
bar_colors = []

for inter in sorted_inters:

    # BOTH
    both = data[inter]['both']
    scatter_points(ax_main, both, current, 'gray',
                   alpha=SCATTER_ALPHA)
    plot_box(ax_main, both, [current], colors['both'])
    nnzero_both = sum(x > 0 for x in both)
    bar_positions.append(current)
    bar_counts.append(nnzero_both)
    bar_colors.append(colors['both'])
    label_pos[inter] = current
    current += 1

    # IMAP only
    imap = data[inter]['imap']
    scatter_points(ax_main, imap, current, 'gray',
                   alpha=SCATTER_ALPHA)
    plot_box(ax_main, imap, [current], colors['imap'])
    nnzero_imap = sum(x > 0 for x in imap)
    bar_positions.append(current)
    bar_counts.append(nnzero_imap)
    bar_colors.append(colors['imap'])
    label_pos[inter] = current
    current += 1

    # PLIF only
    plif = data[inter]['plif']
    scatter_points(ax_main, plif, current, 'gray',
                   alpha=SCATTER_ALPHA)
    plot_box(ax_main, plif, [current], colors['plif'])
    nnzero_plif = sum(x > 0 for x in plif)
    bar_positions.append(current)
    bar_counts.append(nnzero_plif)
    bar_colors.append(colors['plif'])
    label_pos[inter] = current
    current += spacer

bar_positions = np.asarray(bar_positions)
bar_counts = np.asarray(bar_counts)

# Determine where the high axis should start
# ---- set y-lims for broken axis
if len(bar_counts) > 0:
    max_count = bar_counts.max()

    # lower panel: always 0 → LOW_MAX
    ax_top_low.set_ylim(0, LOW_MAX)

    # upper panel: only used if we actually have values above LOW_MAX
    if max_count > LOW_MAX:
        HIGH_MIN = LOW_MAX + 1  # start just above the break
        ax_top_high.set_ylim(HIGH_MIN, max_count * 1.05)
    else:
        # no bar above LOW_MAX; collapse the upper axis
        ax_top_high.set_ylim(0, 1)


# ---- draw bars:
#  - if c <= LOW_MAX: full bar in low axis only
#  - if c > LOW_MAX: truncated bar (LOW_MAX) in low axis + full bar in high axis
for pos, c, col in zip(bar_positions, bar_counts, bar_colors):
    if c <= LOW_MAX:
        ax_top_low.bar(pos, c, color=col, alpha=BAR_ALPHA,
                       width=1, edgecolor='none')
    else:
        # truncated bar in low axis
        ax_top_low.bar(pos, LOW_MAX, color=col, alpha=BAR_ALPHA,
                       width=1, edgecolor='none')
        # full bar in high axis
        ax_top_high.bar(pos, c, color=col, alpha=BAR_ALPHA,
                        width=1, edgecolor='none')



# Diagonal marks for the break
d = 0.015
kwargs = dict(marker=[(-1, -d), (1, d)],
              markersize=8,
              linestyle="none",
              color='k',
              mec='k',
              mew=1,
              clip_on=False)
ax_top_high.plot([0, 1], [0, 0], transform=ax_top_high.transAxes, **kwargs)
ax_top_low.plot([0, 1], [1, 1], transform=ax_top_low.transAxes, **kwargs)

# ---- x-ticks on main panel
pos = np.asarray(list(label_pos.values())) - 1
ax_main.set_xticks(pos)
ax_main.set_xticklabels(
    list(label_pos.keys()),
    rotation=55, ha='right', fontsize=15
)
ax_main.set_xlim(-2, current + 1)

# ---- add centered y-label for the combined top block
# (instead of attaching it only to the lower axis)
bbox_high = ax_top_high.get_position()
bbox_low = ax_top_low.get_position()
y_center_top = 0.5 * (bbox_high.y0 + bbox_low.y1)

fig1.text(0.02, y_center_top, '# Interacting Residue Pairs',
          va='center', rotation='vertical')

# ---- legend anchored to top_high
legend_boxes = [
    plt.Rectangle((0, 0), 1, 1, color=colors['both'], alpha=BAR_ALPHA),
    plt.Rectangle((0, 0), 1, 1, color=colors['imap'], alpha=BAR_ALPHA),
    plt.Rectangle((0, 0), 1, 1, color=colors['plif'], alpha=BAR_ALPHA)
]

ax_top_high.legend(
    legend_boxes,
    ['Both', 'InterMap only', 'ProLIF only'],
    loc='upper center', ncol=3, fontsize='small',
    frameon=False, bbox_to_anchor=(0.5, 1.30)
)

plt.tight_layout()
plt.savefig(join(out_dir, 'identity-per-type.png'),
            bbox_inches='tight')
plt.show()
# plt.close()
# %%

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
# fig1, ax_main = plt.subplots(dpi=300)
# ax_main.set_xlabel('Interaction type', fontweight='bold')
# ax_main.set_ylabel('Interactions detected (%)', fontweight='bold')
# ax_main.set_ylim(0, 105)
#
# for i, inter in enumerate(percentage):
#     both = percentage[order[i][0]]['both'][0]
#     plif = percentage[order[i][0]]['plif'][0]
#     imap = percentage[order[i][0]]['imap'][0]
#     ax_main.bar(i, both, color=cmp.c3, zorder=1, lw=0.5, alpha=cmp.alpha1)
#     ax_main.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1,
#             alpha=cmp.alpha1)
#     ax_main.bar(i, plif, color=cmp.c1, bottom=both, zorder=1, alpha=cmp.alpha1)
# plt.show()
# ax_main.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.5)
# ax_main.bar(i, both, color=cmp.c3, label='Both', alpha=0)
# ax_main.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif,
#         alpha=0)
# ax_main.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both, alpha=0)
#
# leg = ax_main.legend(loc='lower center', ncol=3, fontsize='medium', fancybox=False,
#                  framealpha=0.95)
# for lh in leg.legend_handles:
#     lh.set_alpha(1)
# ax_main.set_xticks(range(len(order)))
# ax_main.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
# plt.tight_layout()
# plt.savefig(join(out_dir, 'identity-per-type-global.png'))
# plt.close()
#
# # %% ==========================================================================
# # Plot 2: Detailed Insights
# # =============================================================================
# fig, ax_main = plt.subplots(1, 1, sharex=True)
#
# # ax2.set_xlabel('Interaction type', fontweight='bold')
# fig.text(0.0, 0.6, '# Detected interactions', va='center', rotation='vertical',
#          fontweight='bold')
#
# for i, inter in enumerate(condensed):
#     both = condensed[order[i][0]]['both'][0]
#     plif = condensed[order[i][0]]['plif'][0]
#     imap = condensed[order[i][0]]['imap'][0]
#     ax_main.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
#     # ax2.bar(i, both, color=cmp.c3, zorder=1, lw=0.5)
#     ax_main.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
#     # ax2.bar(i, imap, color=cmp.c2, bottom=both + plif, zorder=1)
#     ax_main.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)
#     # ax2.bar(i, plif, color=cmp.c1, bottom=both, zorder=1)
#
# # ax_main.set_ylim(0, 100)  # outliers only
# # ax2.set_ylim(0, 10)  # most of the data
#
# ax_main.spines.bottom.set_visible(False)
# # ax2.spines.top.set_visible(False)
#
# ax_main.grid(axis='y', lw=1, ls='--', color='gray', zorder=5, alpha=cmp.alpha3)
# # ax2.grid(axis='y', lw=1, ls='--', color=cmp.c2, zorder=5, alpha=cmp.alpha3)
#
# ax_main.bar(i, both, color=cmp.c3, label='Both')
# # ax2.bar(i, both, color=cmp.c3, label='Both')
# ax_main.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
# # ax2.bar(i, imap, color=cmp.c2, label='InterMap only', bottom=both + plif)
# ax_main.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
# # ax2.bar(i, plif, color=cmp.c1, label='ProLIF only', bottom=both)
# ax_main.legend(loc='upper center', ncol=3, fontsize='medium', fancybox=False,
#            framealpha=cmp.alpha1)
#
# # ax2.set_xticks(range(len(order)))
# # ax2.set_xticklabels([x[0] for x in order], rotation=45, ha='right')
#
# d = 0.25  # proportion of vertical to horizontal extent of the slanted line
# kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
#               linestyle="none", color=cmp.c1, mec='k', mew=1.5, clip_on=False)
# ax_main.plot([0, 1], [0, 0], transform=ax_main.transAxes, **kwargs)
# # ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
#
# plt.tight_layout()
# plt.savefig(join(out_dir, 'identity-per-type-detailed.png'))
# plt.close()
