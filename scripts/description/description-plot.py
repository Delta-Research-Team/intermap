# Created by gonzalezroy at 4/9/25
import os
from os.path import basename, join

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


def plot_img(ax, img_path, title):
    """
    Plot an image in a given axis with a title.

    Args:
        ax: axis object
        img_path: path to the image
        title: title of the image
    """
    ax.imshow(plt.imread(img_path))
    ax.set_title(title, fontsize=16, fontname='STIXGeneral')
    ax.axis('off')


# =============================================================================
#
# =============================================================================
user = os.getenv('USER')
img_dir = f'/home/{user}/RoyHub/intermap/scripts/description'
images = {basename(x): x for x in gnl.recursive_finder('*.png', img_dir)}
out_dir = f'/home/{user}/RoyHub/intermap/scripts/description'

# =============================================================================
#
# =============================================================================
os.makedirs(out_dir, exist_ok=True)
import matplotlib.colors as colors

df = pd.DataFrame.from_dict(systems)
bounds = np.arange(4, df.max().max(), 10000)
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)

# >>>> Layout
fig = plt.figure(constrained_layout=False, figsize=(14, 8), dpi=600)
gs = fig.add_gridspec(2, 5, hspace=0.2, wspace=0, width_ratios=[1, 1, 1, 1, 1],
                      height_ratios=[1, 1])

mpro = fig.add_subplot(gs[0, 0])
igg = fig.add_subplot(gs[0, 1], sharey=mpro)
p53 = fig.add_subplot(gs[0, 2], sharey=mpro)
nsp13 = fig.add_subplot(gs[0, 3], sharey=mpro)
pao1 = fig.add_subplot(gs[0, 4])
spike = fig.add_subplot(gs[1, 0], sharey=pao1)
ace2 = fig.add_subplot(gs[1, 1], sharey=pao1)
mdlr = fig.add_subplot(gs[1, 2], sharey=mpro)
nucleosome = fig.add_subplot(gs[1, 3], sharey=pao1)

plot_img(mpro, images['mpro1.png'], 'MPRO')
plot_img(igg, images['M1_IgG3.png'], 'IgG3-M1')
plot_img(p53, images['p53.png'], 'p53')
plot_img(nsp13, images['nsp13.png'], 'NSP13')
plot_img(pao1, images['PAO1_prot.png'], 'EC_T4P')
plot_img(spike, images['spike.png'], 'Spike')
plot_img(ace2, images['ACE2.png'], 'ACE2-RBD')
plot_img(mdlr, images['mdlr.png'], 'Topoiso_I-CPT')
plot_img(nucleosome, images['nucleosome.png'], 'Nucleosome')
fig.align_titles()
# plt.subplots_adjust(left=0, right=1, top=1, bottom=0, hspace=0, wspace=0)
fig.savefig(join(out_dir, 'systems.png'), dpi=300, bbox_inches='tight')
plt.close()
