# Created by rglez at 4/6/25
import re
from os.path import join

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rgpack import generals as gnl

import intermap.commonplots as cmp


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


def plot_block(matrix, axis, cmap, vmin=None, vmax=None):
    im1 = axis.imshow(matrix, cmap=cmap, aspect='auto',
                      vmin=vmin, vmax=vmax,
                      alpha=cmp.alpha1)

    axis.set_yticks(range(len(matrix.index)))
    axis.set_yticklabels(matrix.index)
    axis.set_yticks(np.arange(-.5, matrix.shape[0], 1), minor=True)
    axis.set_xticks(range(len(matrix.columns)))
    axis.set_xticklabels(matrix.columns)
    axis.set_xticks(np.arange(-.5, matrix.shape[1], 1), minor=True)
    axis.tick_params(axis='y', which='minor', length=0)
    axis.tick_params(axis='x', which='minor', length=0)
    axis.tick_params(axis='x', which='major', length=0)
    axis.grid(True, axis='both', color='k', linestyle='-',
              linewidth=cmp.lw1, which='minor')
    return im1


def plot_colorbar(ax1, ax2, im1, title=None):
    cax, kw = matplotlib.colorbar.make_axes([ax1, ax2])
    cb = plt.colorbar(im1, cax=cax, **kw)
    for l in cb.ax.yaxis.get_ticklabels():
        l.set_fontsize(cmp.fs1)
        l.set_color('k')
    cb.ax.set_ylabel(title, fontsize=cmp.fs1, labelpad=2, fontweight='regular')
    cb.ax.minorticks_on()


# =============================================================================
# User-defined variables
# =============================================================================
user = 'gonzalezroy'
root_dir = '/home/rglez/RoyHub/intermap/scripts/scalability'
out_dir = root_dir
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


# >>>> Grid layout
gridspec = dict(hspace=0., wspace=0., height_ratios=[1, 0.1, 1])
fig, axs = plt.subplots(
    nrows=3, ncols=2, gridspec_kw=gridspec, sharex=True, sharey=True,
    constrained_layout=True, figsize=(8, 8), dpi=300)

ax_atom_time = axs[0, 1]
ax_resid_time = axs[0, 0]
axs[1, 0].set_visible(False)
axs[1, 1].set_visible(False)
ax_atom_ram = axs[2, 1]
ax_resid_ram = axs[2, 0]

ax_atom_time.set_title('Resolution: "atom"', fontweight='bold')
ax_atom_ram.set_title('Resolution: "atom"', fontweight='bold')
ax_resid_ram.set_title('Resolution: "residue"', fontweight='bold')
ax_resid_time.set_title('Resolution: "residue"', fontweight='bold')

# >>>> Prepare data
atom_time = atomic.pivot(index='chunk_size', columns='n_procs',
                         values='time')
resid_time = residic.pivot(index='chunk_size', columns='n_procs',
                           values='time')
atom_ram = atomic.pivot(index='chunk_size', columns='n_procs',
                        values='ram')
resid_ram = residic.pivot(index='chunk_size', columns='n_procs',
                          values='ram')

# Plot time
minim = min(atom_time.min().min(), resid_time.min().min())
maxim = max(atom_time.max().max(), resid_time.max().max())
im1 = plot_block(atom_time, ax_atom_time, cmap=cmp.cmap1, vmin=minim,
                 vmax=maxim)
im2 = plot_block(resid_time, ax_resid_time, cmap=cmp.cmap1, vmin=minim,
                 vmax=maxim)
plot_colorbar(ax_atom_time, ax_resid_time, im1, title='Time (min)')

# Plot ram
minim = min(atom_ram.min().min(), resid_ram.min().min())
maxim = max(atom_ram.max().max(), resid_ram.max().max())
im1 = plot_block(atom_ram, ax_atom_ram, cmap=cmp.cmap1, vmin=minim, vmax=maxim)
im2 = plot_block(resid_ram, ax_resid_ram, cmap=cmp.cmap1, vmin=minim,
                 vmax=maxim)

plot_colorbar(ax_atom_ram, ax_resid_ram, im1, title='RAM (GB)')

# Set labels
ax_atom_ram.set_xlabel('# Processors', fontweight='roman')
ax_resid_ram.set_xlabel('# Processors', fontweight='roman')

ax_resid_ram.set_ylabel('Chunk size', fontweight='roman')
ax_resid_time.set_ylabel('Chunk size', fontweight='roman')

plt.savefig(join(out_dir, 'scalability.png'))
plt.close()
