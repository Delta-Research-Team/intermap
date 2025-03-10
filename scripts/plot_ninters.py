# Created by rglez at 3/8/25
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from rgpack import generals as gnl


def get_data(data, traj, sele):
    try:
        return data[traj][sele]
    except KeyError:
        return 0


def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))


# =============================================================================
# Get raw data
data_plif = '/home/rglez/RoyHub/fault/intermap/scripts/prolif_data.pickle'
data_gc = '/home/rglez/RoyHub/fault/intermap/scripts/gc_data.pickle'
data_imap = '/home/rglez/RoyHub/fault/intermap/scripts/imap_data.pickle'
out_dir = '/home/rglez/RoyHub/fault/intermap/scripts/'

imap = gnl.unpickle_from_file(data_imap)
gc = gnl.unpickle_from_file(data_gc)
plif = gnl.unpickle_from_file(data_plif)

# Get trajectories and selections
trajs = list(imap.keys())
seles = set()
for traj in imap:
    for sele in imap[traj]:
        seles.add(sele)

# Get data  for each software
imap_data = defaultdict(lambda: defaultdict(int))
gc_data = defaultdict(lambda: defaultdict(int))
plif_data = defaultdict(lambda: defaultdict(int))
for traj in trajs:
    for sele in seles:
        pass
        try:
            imap_data[traj][sele] = sum(imap[traj][sele].values())
        except KeyError:
            imap_data[traj][sele] = 0

        try:
            gc_data[traj][sele] = sum(gc[traj][sele].values())
        except KeyError:
            gc_data[traj][sele] = 0

        try:
            plif_data[traj][sele] = sum(plif[traj][sele].values())
        except KeyError:
            plif_data[traj][sele] = 0

font = {'family': 'sans-serif',
        'size': 12}
matplotlib.rc('font', **font)

traj_names = {'1kx5sno_dry': 'Nucl_SNO',
              '8oxoGA2_1': 'Nucl_8OXO',
              'M1': 'IgG-M1',
              'ncp-OCS-nosalt': 'Nucl_OCS',
              'prolif_tutorial': 'GPCR-Lig'}

barWidth = 0.25
br1 = np.arange(len(trajs))
br2 = [x + barWidth for x in br1]
br3 = [x + barWidth for x in br2]

for sele in seles:
    IM = [get_data(imap_data, x, sele) for x in trajs]
    IM_percent = [(x / x) * 100 if x > 0 else 0 for x in IM]

    GC = [get_data(gc_data, x, sele) for x in trajs]
    GC_percent = [round((GC[i] / IM[i]) * 100, 2) if x > 0 else 0 for i, x in
                  enumerate(IM)]

    PL = [get_data(plif_data, x, sele) for x in trajs]
    PL_percent = [round((PL[i] / IM[i]) * 100, 2) if x > 0 else 0 for i, x in
                  enumerate(IM)]

    # Starting the plot
    fig, ax1 = plt.subplots(dpi=600)
    fig.suptitle(sele)
    ax1.set_xlabel('System label', fontweight='bold')
    ax1.set_ylabel('% of detected interactions', fontweight='bold')
    ax1.set_ylim(0, 120)

    ax1.bar(br1, GC_percent, color='k', width=barWidth, label='getContacts')
    ax1.bar(br2, PL_percent, color='darkgrey', width=barWidth, label='ProLIF')
    ax1.bar(br3, IM_percent, color='grey', width=barWidth, label='InterMap')

    for i, txt in enumerate(IM):
        if IM_percent[i] > 0:
            ax1.text(br3[i], 103,
                     r"${0:s}$".format(as_si(txt, 1)), ha='center',
                     va='bottom', fontweight='bold', fontsize=10, color='grey')
    ax1.set_xticks(np.arange(len(trajs)) + barWidth)
    ax1.set_xticklabels([traj_names[x] for x in trajs], rotation=45, ha='right')
    ax1.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.25)

    ax1.legend(loc='center', ncol=3, fontsize='medium', fancybox=False,
               framealpha=0.2, bbox_to_anchor=(0.5, 1.1))
    plt.tight_layout()


    plt.savefig(out_dir + sele + '.png')
    plt.close()
    print(sele)
