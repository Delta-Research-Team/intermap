# Created by gonzalezroy at 3/26/25
import matplotlib.pyplot as plt
import numpy as np

import benchmark_data as bd

# =============================================================================
#
# =============================================================================
info = bd.info
n_seles = 0
n_trajs = 0
for s, software in enumerate(info):
    for traj in info[software]:
        n_trajs += 1
        for sele in info[software][traj]:
            n_seles += 1

num_soft = len(info)
num_sele = n_seles // num_soft
num_traj = n_trajs // num_soft

ram_array = np.zeros((num_sele, num_soft), dtype=float)
time_array = ram_array.copy()
n_inters_array = ram_array.copy()

for s, software in enumerate(info):
    line = -1
    sele_label = []
    sele_limits = []
    for t, traj in enumerate(info[software]):
        for l, sele in enumerate(info[software][traj]):
            line += 1
            sele_label.append(sele)
            sele_limits.append(traj)
            ram = info[software][traj][sele]['ram']
            time = info[software][traj][sele]['time']
            n_inters = info[software][traj][sele]['n_inters']
            ram_array[line, s] = ram
            time_array[line, s] = time
            n_inters_array[line, s] = np.random.randint(1, 100000000)

limits_raw = np.asarray(sele_limits)
limits_int = np.arange(0, len(limits_raw))
limits = []
traj_labels = []
for traj in np.unique(limits_raw):
    indices = np.where(limits_raw == traj)[0]
    traj_labels.append(traj)
    limits.append(limits_int[indices].mean())

# =============================================================================
#
# =============================================================================
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_block(array, axis, cmap, cbar_label):
    im1 = axis.matshow(array, cmap=cmap, norm=LogNorm())
    axis.set_xticks(range(num_soft), labels=list(info.keys()), rotation=45,
                    ha='left')
    axis.set_yticks(range(num_sele), labels=sele_label)

    axis.hlines(limits, -0.5, num_soft - 0.5, color='k', lw=1)

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='10%', pad=0.0)
    cb = fig.colorbar(im1, cax=cax, orientation='vertical', label='Time (s)')
    cb.ax.set_ylabel(cbar_label, font="Ubuntu mono", fontsize=16,
                     labelpad=5, color='k', fontweight='bold')

    for l in cb.ax.yaxis.get_ticklabels():
        l.set_family("Ubuntu mono")
        l.set_fontsize(14)
        l.set_color('k')


cmap = 'cividis'
fig, (ax_time, ax_ram, ax_nints) = plt.subplots(1, 3, figsize=(12, 6),
                                                dpi=600, )
plot_block(time_array, ax_time, cmap, 'Wallclock Time (s)')
plot_block(ram_array, ax_ram, cmap, 'RAM Peak (GB)')
plot_block(n_inters_array, ax_nints, cmap, '# Detected Interactions')
plt.tight_layout()
plt.subplots_adjust(right=0.99)
plt.savefig('benchmark_plot.png')
plt.close()
