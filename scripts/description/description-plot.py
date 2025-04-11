# Created by gonzalezroy at 4/9/25
from os.path import basename

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rgpack import generals as gnl

systems = {
    'MPRO': {
        'Protein': 4674,
        'DNA': 0,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 0,
        'Metal': 0,
        'Total': 4674,
        "Sel-1": 4674,
        "Sel-2": 4674,
        'Frames': 110000,
    },

    'IgG3-M1': {
        "Protein": 26502,
        'DNA': 0,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 0,
        'Metal': 0,
        'Total': 26502,
        "Sel-1": 4718,
        "Sel-2": 21784,
        'Frames': 30000,
    },

    'p53': {
        "Protein": 24178,
        'DNA': 3177,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 0,
        'Metal': 4,
        'Total': 27329,
        "Sel-1": 3177,
        "Sel-2": 24148,
        'Frames': 99000,
    },

    'NSP13': {
        "Protein": 40027,
        'DNA': 2863,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 0,
        'Metal': 8,
        'Total': 42898,
        "Sel-1": 2871,
        "Sel-2": 40027,
        'Frames': 1000,
    },

    'PAO1': {
        "Protein": 46376,
        'DNA': 0,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 0,
        'Metal': 0,
        'Total': 46376,
        "Sel-1": 46376,
        "Sel-2": 46376,
        'Frames': 20000,
    },

    'Spike': {
        "Protein": 58523,
        'DNA': 0,
        'Glycan': 14236,
        'Lipid': 0,
        'Water': 0,
        'Metal': 0,
        'Total': 72759,
        "Sel-1": 14236,
        "Sel-2": 58523,
        'Frames': 45000,
    },

    'ACE2-RBD': {
        'Protein': 30884,
        'DNA': 0,
        'Glycan': 4068,
        'Lipid': 109311,
        'Water': 0,
        'Metal': 0,
        'Total': 144263,
        "Sel-1": 30884,
        "Sel-2": 113379,
        'Frames': 30000,
    },

    'Nucleosome': {
        'Protein': 15729,
        'DNA': 9346,
        'Glycan': 0,
        'Lipid': 0,
        'Water': 372543,
        'Metal': 0,
        'Total': 397776,
        "Sel-1": 9346,
        "Sel-2": 388272,
        'Frames': 20757,
    }
}

# =============================================================================
#
# =============================================================================
img_dir = '/home/gonzalezroy/RoyHub/intermap/data/figures/systems/'
images = {basename(x): x for x in gnl.recursive_finder('*.png', img_dir)}

# =============================================================================
#
# =============================================================================
import matplotlib.colors as colors

df = pd.DataFrame.from_dict(systems)
bounds = np.arange(4, df.max().max(), 10000)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# >>>> Layout
fig = plt.figure(constrained_layout=True, figsize=(12, 8), dpi=600)
gs = fig.add_gridspec(2, 4, hspace=0., wspace=0, width_ratios=[1, 1, 1, 1],
                      height_ratios=[1, 1])


def plot_img(ax, img_path, title):
    ax.imshow(plt.imread(img_path))
    ax.set_title(title, fontweight='bold', fontsize=20, fontname='Ubuntu mono')
    ax.axis('off')
    pass


plot_img(fig.add_subplot(gs[0, 0]), images['mpro1.png'], 'MPRO')
plot_img(fig.add_subplot(gs[0, 1]), images['M1_IgG3.png'], 'IgG3-M1')
plot_img(fig.add_subplot(gs[0, 2]), images['p53.png'], 'p53')
plot_img(fig.add_subplot(gs[0, 3]), images['nsp13.png'], 'NSP13')
plot_img(fig.add_subplot(gs[1, 0]), images['PAO1_prot.png'], 'PAO1')
plot_img(fig.add_subplot(gs[1, 1]), images['spike.png'], 'Spike')
plot_img(fig.add_subplot(gs[1, 2]), images['ACE2.png'], 'ACE2-RBD')
plot_img(fig.add_subplot(gs[1, 3]), images['nucleosome.png'], 'Nucleosome')

# inv = fig.add_subplot(gs[2, :])
# inv.axis('off')

# ax = fig.add_subplot(gs[3, :])
# ax2 = fig.add_subplot(gs[4, :])
# ax.invert_yaxis()
# ax.xaxis.tick_top()

# im = ax.imshow(df.T, cmap='cividis', aspect='auto', alpha=0.95, norm=norm)
#
# ax.set_yticks(range(len(df.columns)), labels=df.columns,
#               fontname='Ubuntu mono', fontsize=12)
# ax.set_yticks(np.arange(-.5, len(df.columns), 1), minor=True)
# ax.set_xticks(range(len(df.index)), labels=df.index, rotation=0, ha='center',
#               fontname='Ubuntu mono', fontsize=12)
# ax.set_xticks(np.arange(-.5, len(df.index), 1), minor=True)
# ax.tick_params(axis='both', which='minor', length=0)
# ax.grid(color='k', lw=0.5, alpha=0.5, which='minor', axis='both', ls='-')
#
# ax.axvline(x=6.5, color='k', lw=1)
# cb = fig.colorbar(im, cax=ax2, orientation='horizontal')
# cb.ax.xaxis.set_ticks_position('bottom')
#
# cb.ax.set_xlabel('[Atoms | Frames]', font="Ubuntu mono", fontsize=12,
#                  labelpad=5, color='k', fontweight='bold')
#
# for l in cb.ax.yaxis.get_ticklabels():
#     l.set_family("Ubuntu mono")
#     l.set_fontsize(8)
#     l.set_color('k')

plt.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)
fig.savefig('systems.png', dpi=300, bbox_inches='tight')
plt.close()
