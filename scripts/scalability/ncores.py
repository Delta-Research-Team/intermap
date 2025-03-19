# Created by gonzalezroy at 3/11/25
from os.path import basename

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import rgpack.generals as gnl


def get_ram_time(file):
    with open(file) as f:
        lines = f.readlines()
    ram, time = 0, 0
    for line in lines:
        if line.startswith('\tUser time'):
            try:
                time = float(line.split()[3])
            except ValueError:
                time = 0
        if line.startswith('\tMaximum resident set size'):
            try:
                ram = int(line.split()[5]) / 1024
            except ValueError:
                ram = 0
    return ram, time


def get_perfs(files, sep):
    perfs = gnl.recursive_defaultdict()
    for file in files:
        ncores = int(basename(file).split(sep)[-1])
        ram, time = get_ram_time(file)
        perfs[ncores]['ram'] = ram
        perfs[ncores]['time'] = time
    return perfs


def plot_perfs(perfs, range2plot, out_name, xlabel):
    font = {'family': 'sans-serif', 'size': 12}
    matplotlib.rc('font', **font)
    barWidth = 0.25
    br1 = np.arange(len(perfs))
    br2 = [x + barWidth for x in br1]

    # Starting the plot
    fig, ax1 = plt.subplots(dpi=600, tight_layout=True)
    ax2 = ax1.twinx()
    ax1.set_xlabel('No. of cores', fontweight='bold')
    ax1.set_ylabel('Time (s)', fontweight='bold')
    ax2.set_ylabel('RAM (MB)', fontweight='bold', color='darkgrey')
    ax2.spines['right'].set_color('darkgrey')

    # Create bars and assign labels directly
    bar1 = ax1.bar(br1, [perfs[x]['time'] for x in range2plot], color='k',
                   width=barWidth, label='Time')
    bar2 = ax2.bar(br2, [perfs[x]['ram'] for x in range2plot], color='darkgrey',
                   width=barWidth, label='RAM')

    ax1.set_xticks(np.arange(len(perfs)))
    ax1.set_xticklabels(list(range2plot), ha='right', rotation=45)
    ax1.grid(axis='y', lw=1, ls='--', color='k', zorder=5, alpha=0.25)

    ax2.tick_params(axis='y', colors='darkgrey')

    # Create legends for both axes
    ax1.legend(handles=[bar1], loc='upper left', fontsize='medium', fancybox=False,
               framealpha=0, bbox_to_anchor=(0.5, 1.1))
    ax2.legend(handles=[bar2], loc='upper right', fontsize='medium', fancybox=False,
               framealpha=0, bbox_to_anchor=(0.5, 1.1))

    plt.savefig(out_name)
    plt.close()


# =============================================================================
#
# =============================================================================
path_cores = '/media/gonzalezroy/Expansion/RoyData/intermap/M1/core/'
path_chunks = '/media/gonzalezroy/Expansion/RoyData/intermap/M1/chunk/'
out_dir = '/home/gonzalezroy/RoyHub/intermap/scripts/'

files_cores = list(gnl.recursive_finder('core*', path_cores))
perfs_cores = get_perfs(files_cores, sep='core')
range_cores = sorted(perfs_cores.keys())
out_name = out_dir + 'ncores.png'
plot_perfs(perfs_cores, range_cores, out_name, 'No. of cores')

files_chunks = list(gnl.recursive_finder('chunk*', path_chunks))
perfs_chunks = get_perfs(files_chunks, sep='chunk')
range_chunks = sorted(perfs_chunks.keys())
out_name = out_dir + 'nchunks.png'
plot_perfs(perfs_chunks, range_chunks, out_name, 'No. of chunks')

# =============================================================================
