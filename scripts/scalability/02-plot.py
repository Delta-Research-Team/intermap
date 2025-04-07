# Created by rglez at 4/6/25
import re

import pandas as pd
from rgpack import generals as gnl


def parse_performance(txt):
    # Regular expression patterns
    elapsed_time = r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\d+:\d+\.\d+)'
    max_memory = r'Maximum resident set size \(kbytes\): (\d+)'
    exit_status = r'Command terminated by signal (\d+)'

    # Read the file
    with open(txt) as f:
        text = f.read()

    # Finding matches
    elapsed_time_match = re.findall(elapsed_time, text)[0]
    max_memory_match = re.findall(max_memory, text)[0]
    exit_status_match = re.findall(exit_status, text)

    # Process matches
    elapsed_time = elapsed_time_match.split(':')
    if len(elapsed_time) == 2:
        time = int(elapsed_time[0]) * 60 + float(elapsed_time[1])
    else:
        time = int(elapsed_time[0]) * 3600 + int(elapsed_time[1]) * 60 + float(
            elapsed_time[2])

    ram = int(max_memory_match) / 2 ** 20
    status = int(exit_status_match[0]) if exit_status_match else 0
    return round(time / 60, 2), round(ram, 2), status


# =============================================================================
# User-defined variables
# =============================================================================
root_dir = '/home/rglez/RoyHub/intermap/data/scalability_by_RGA'
# =============================================================================

# Find all .txt files in the root directory
txts = list(gnl.recursive_finder('*.txt', root_dir))

# Parse the performance data
data = []
for txt in txts:
    res, proc, chunk = txt.split('/')[-2].split('-')
    time, ram, status = parse_performance(txt)
    data.append((res, int(proc), int(chunk), time, ram, status))

df = pd.DataFrame(data, columns=['resolution', 'n_procs', 'chunk_size', 'time',
                                 'ram', 'status'])
df.sort_values(by=['resolution', 'n_procs', 'chunk_size'], inplace=True)

# Rearrange
to_plot = df[df['status'] == 0]
to_plot.reset_index(drop=True, inplace=True)
atomic = to_plot[to_plot['resolution'] == 'atom']
residic = to_plot[to_plot['resolution'] == 'residue']

# %% ==========================================================================
# Plotting
# =============================================================================
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Plot constants
font = {'family': 'sans-serif',
        'size': 12}
matplotlib.rc('font', **font)

# Grid layout
gridspec = dict(hspace=0.15, height_ratios=[1, 1, 0.25, 0.25, 1, 1])
fig, axs = plt.subplots(nrows=3, ncols=2, gridspec_kw=gridspec, sharex=True,
                        figsize=(4, 8), dpi=300)

axs[2].set_visible(False)
axs[3].set_visible(False)

plt.show()

atom_time = axs[0]
atom_ram = axs[3]
resid_time = axs[1]
resid_ram = axs[4]

# Plot atomic data
atom_matrix_time = atomic.pivot(index='chunk_size', columns='n_procs',
                                values='time')
atom_matrix_ram = atomic.pivot(index='chunk_size', columns='n_procs',
                               values='ram')
resid_matrix_time = residic.pivot(index='chunk_size', columns='n_procs',
                                  values='time')
resid_matrix_ram = residic.pivot(index='chunk_size', columns='n_procs',
                                 values='ram')


def plot_block(matrix, axis, cmap):
    im1 = axis.imshow(matrix, cmap=cmap, aspect='auto')
    axis.set_yticks(range(len(matrix.index)))
    axis.set_yticklabels(matrix.index)
    axis.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    axis.set_xticks(range(len(matrix.columns)))
    axis.set_xticklabels(matrix.columns)
    axis.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    axis.tick_params(axis='y', which='minor', length=0)
    axis.tick_params(axis='x', which='minor', length=0)
    axis.tick_params(axis='x', which='major', length=0)
    axis.grid(True, axis='both', color='white', linestyle='--', linewidth=0.5,
              which='minor')
    return im1


cmap = 'cividis'
im1 = plot_block(atom_matrix_time, atom_time, cmap=cmap)
im2 = plot_block(resid_matrix_time, resid_time, cmap=cmap)

cax, kw = matplotlib.colorbar.make_axes([atom_time, resid_time], )
cb = plt.colorbar(im1, cax=cax, **kw)
for l in cb.ax.yaxis.get_ticklabels():
    l.set_family("Ubuntu mono")
    l.set_fontsize(14)
    l.set_color('k')
cb.ax.set_ylabel('Time (min)', font="Ubuntu mono", fontsize=16,
                 labelpad=5, color='k', fontweight='bold')

im3 = plot_block(atom_matrix_ram, atom_ram, cmap=cmap)
im4 = plot_block(resid_matrix_ram, resid_ram, cmap=cmap)

cax2, kw2 = matplotlib.colorbar.make_axes([atom_ram, resid_ram], )
cb2 = plt.colorbar(im3, cax=cax2, **kw2)
for l in cb2.ax.yaxis.get_ticklabels():
    l.set_family("Ubuntu mono")
    l.set_fontsize(14)
    l.set_color('k')
cb2.ax.set_ylabel('RAM (GB)', font="Ubuntu mono", fontsize=16,
                  labelpad=5, color='k', fontweight='bold')

plt.tight_layout()
plt.savefig('scalability.png')
plt.close()
plt.show()
