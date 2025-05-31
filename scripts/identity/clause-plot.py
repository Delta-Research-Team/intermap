# Created by rglez at 5/24/25

from os.path import join

import matplotlib.pyplot as plt


def create_nature_quality_boxplot(data, out_dir):
    """
    Create a publication-ready boxplot for Nature journal standards
    """
    # Set publication-ready style
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 12,
        'axes.linewidth': 1.5,
        'xtick.major.width': 1.5,
        'ytick.major.width': 1.5,
        'xtick.minor.width': 1,
        'ytick.minor.width': 1,
        'legend.frameon': False,
        'axes.spines.top': False,
        'axes.spines.right': False
    })

    n_inters = len(data)

    # Create figure with appropriate size and DPI for publication
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

    # Define colors (using colorblind-friendly palette)
    colors = {
        'both': '#2E86AB',  # Blue
        'imap': '#A23B72',  # Magenta
        'plif': '#F18F01'  # Orange
    }

    # Create positions for grouped boxplots
    positions = []
    labels = []
    x_labels = []

    group_width = 3
    current_pos = 0

    all_data = []
    all_positions = []
    all_colors = []

    for i, inter in enumerate(data.keys()):
        # Collect data for each method
        methods = ['both', 'imap', 'plif']
        group_positions = [current_pos + j for j in range(len(methods))]

        for j, method in enumerate(methods):
            all_data.append(data[inter][method])
            all_positions.append(group_positions[j])
            all_colors.append(colors[method])

        # Store center position for x-axis labels
        x_labels.append(inter)
        labels.append(current_pos + 1)  # Center of the group

        current_pos += group_width + 1

    # Create boxplots
    bp = ax.boxplot(all_data, positions=all_positions, patch_artist=True,
                    widths=0.6, showfliers=True,
                    boxprops=dict(linewidth=1.5),
                    whiskerprops=dict(linewidth=1.5),
                    capprops=dict(linewidth=1.5),
                    medianprops=dict(linewidth=2, color='white'),
                    flierprops=dict(marker='o', markersize=4, alpha=0.6))

    # Color the boxes
    for patch, color in zip(bp['boxes'], all_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        patch.set_edgecolor('black')

    # Customize axes
    ax.set_xlabel('Interaction Type', fontsize=14, fontweight='bold',
                  labelpad=10)
    ax.set_ylabel('Interaction Count', fontsize=14, fontweight='bold',
                  labelpad=10)

    # Set x-axis ticks and labels
    ax.set_xticks(labels)
    ax.set_xticklabels(x_labels, fontsize=12)

    # Improve y-axis
    ax.tick_params(axis='y', labelsize=11)
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)

    # Create legend
    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['both'], alpha=0.8,
                      label='Both'),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['imap'], alpha=0.8,
                      label='IMAP'),
        plt.Rectangle((0, 0), 1, 1, facecolor=colors['plif'], alpha=0.8,
                      label='PLIF')]

    ax.legend(handles=legend_elements, loc='upper right', fontsize=11,
              frameon=True, fancybox=True, shadow=True, framealpha=0.9)

    # Add statistical annotations (if needed)
    # You can add significance bars here if you have statistical tests

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Save with high quality for publication
    output_path = join(out_dir, 'identity-per-type.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')

    # Also save as PDF for vector graphics (preferred by many journals)
    pdf_path = join(out_dir, 'identity-per-type.pdf')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')

    # Save as SVG for perfect scalability
    svg_path = join(out_dir, 'identity-per-type.svg')
    plt.savefig(svg_path, format='svg', bbox_inches='tight',
                facecolor='white', edgecolor='none')

    plt.show()

    return fig, ax

# Usage example:
# Assuming your data structure is:
# data = {
#     'interaction1': {'both': [values], 'imap': [values], 'plif': [values]},
#     'interaction2': {'both': [values], 'imap': [values], 'plif': [values]},
#     ...
# }

# Call the function
# fig, ax = create_nature_quality_boxplot(data, out_dir)
