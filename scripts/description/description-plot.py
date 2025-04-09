# Created by gonzalezroy at 4/9/25

systems = {
    'ACE2-RBD': {
        "# Atoms (Sel-1)": 30884,
        "# Atoms (Sel-2)": 113379,
        # '# Total Atoms': 144263,
        '# Frames': 30000,
    },
    'Spike-Open': {
        "# Atoms (Sel-1)": 14236,
        "# Atoms (Sel-2)": 58523,
        # '# Total Atoms': 72759,
        '# Frames': 45000,
    },
    'p53': {
        "# Atoms (Sel-1)": 3177,
        "# Atoms (Sel-2)": 24148,
        # '# Total Atoms': 27329,
        '# Frames': 99000,
    },
    'IgG3-M1': {
        "# Atoms (Sel-1)": 21784,
        "# Atoms (Sel-2)": 4718,
        # '# Total Atoms': 26502,
        '# Frames': 30000,
    },
    'PAO1': {
        "# Atoms (Sel-1)": 46376,
        "# Atoms (Sel-2)": 46376,
        # '# Total Atoms': 46376,
        '# Frames': 20000,
    },
    'Nucleosome': {
        "# Atoms (Sel-1)": 9346,
        "# Atoms (Sel-2)": 15729,
        # '# Total Atoms': 397776,
        '# Frames': 20757,
    },
    'MPRO': {
        "# Atoms (Sel-1)": 4674,
        "# Atoms (Sel-2)": 4674,
        # '# Total Atoms': 4674,
        '# Frames': 110000,
    },
}

# =============================================================================
#
# =============================================================================
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import matplotlib.ticker

# Set the layout
n_cases = len(systems)
n_info = 4
n_spacing = 4
x_len = (n_cases * n_info) + n_spacing * (n_cases - 1)
color_info = {'# Atoms (Sel-1)': 'blue',
              '# Atoms (Sel-2)': 'orange',
              '# Total Atoms': 'green',
              '# Frames': 'purple'}

# Create a list of x_labels, x_values, and x_names
x_labels = []
x_values = []
x_names = []
for x in systems.keys():
    for y in systems[x].keys():
        x_labels.append(y)
        x_values.append(systems[x][y])
        x_names.append(x)
    for z in range(n_spacing):
        x_labels.append('')
        x_values.append(0)
        x_names.append('')
x_labels = np.asarray(x_labels)
x_values = np.asarray(x_values)

# Create a color array based on the x_labels
colors = [color_info[x] if x in color_info.keys() else 'white'
          for x in x_labels]
colors = np.asarray(colors)

# Separate the indices for the twin axis
whole_idx = range(len(colors))
subset_idx = np.where(x_labels != '# Total Atoms')[0]
twin_idx = np.where(x_labels == '# Total Atoms')[0]

# Create a figure and axis
fig, ax = plt.subplots(figsize=(3, 9), dpi=300)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax.barh(subset_idx, x_values[subset_idx], color=colors[subset_idx],
        alpha=0.7, height=1)
# ax2 = ax.twinx()
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 1))
ax.barh(twin_idx, x_values[twin_idx], color=colors[twin_idx],
        alpha=0.7, height=1)

# set the y-axis equivalent
# l = ax.get_ylim()
# l2 = ax2.get_ylim()
# f = lambda x: l2[0] + (x - l[0]) / (l[1] - l[0]) * (l2[1] - l2[0])
# ticks = f(ax.get_yticks())
# ax2.yaxis.set_major_locator(mpl.ticker.FixedLocator(ticks))

# find unique xticks
x_names = np.asarray(x_names)
x_ticks = []
x_labels = []
for name in set(x_names):
    if name == '':
        continue
    idx = np.where(x_names == name)[0].mean()
    x_ticks.append(idx)
    x_labels.append(name)

x_ticks = np.asarray(x_ticks)
x_labels = np.asarray(x_labels)
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, rotation=0, ha='center')

ax.set_ylabel('# Atoms  |  # Frames', fontweight='bold')
# ax2.set_ylabel('# Total Atoms', color='green', fontweight='bold')
# ax2.tick_params(axis='y', labelcolor='green')

ax.grid(axis='y', linestyle='--', alpha=0.7)
legend_labels = list(color_info.keys())
legend_colors = list(color_info.values())
handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in legend_colors]
ax.legend(handles, legend_labels, loc='upper center', ncol=2)
plt.tight_layout()

# Save the figure
fig.savefig('systems_atoms.png', dpi=300, bbox_inches='tight')
plt.close()
