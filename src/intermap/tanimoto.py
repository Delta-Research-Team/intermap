import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

FILE_PATH = "/home/fajardo/03_Fajardo_Hub/02_InterMap/02_tests/interacciones_prolif.csv"


def calculate_tanimoto_similarity(vector1, vector2):
    intersection = np.sum((vector1 == 1) & (vector2 == 1))
    union = np.sum((vector1 == 1) | (vector2 == 1))
    return intersection / union if union > 0 else 0


def create_similarity_matrix(data, step=10):
    selected_columns = list(range(0, data.shape[1], step))
    n_selected = len(selected_columns)

    similarity_matrix = np.zeros((n_selected, n_selected))

    for i, col_i in enumerate(selected_columns):
        for j, col_j in enumerate(selected_columns):
            vector1 = data.iloc[:, col_i].values
            vector2 = data.iloc[:, col_j].values
            similarity_matrix[i, j] = calculate_tanimoto_similarity(vector1, vector2)

    frame_numbers = [int(data.columns[i]) for i in selected_columns]
    return pd.DataFrame(similarity_matrix, index=frame_numbers, columns=frame_numbers)


def create_similarity_vector(data, step=10):
    selected_columns = list(range(0, data.shape[1], step))

    similarity_vector = []
    frame_numbers = []

    reference_vector = data.iloc[:, 0].values

    for col_idx in selected_columns:
        current_vector = data.iloc[:, col_idx].values
        similarity = calculate_tanimoto_similarity(reference_vector, current_vector)
        similarity_vector.append(similarity)
        frame_numbers.append(int(data.columns[col_idx]))

    return pd.Series(similarity_vector, index=frame_numbers)


def plot_tanimoto_comparison(similarity_series, output_file="tanimoto_comparison.png"):

    plt.figure(figsize=(10, 6), dpi=200)

    plt.plot(similarity_series.index, similarity_series.values,
             marker='o', linestyle='-', linewidth=1.5, markersize=4)

    plt.title("Tanimoto Similarity vs Reference Frame")
    plt.xlabel("Frame Number")
    plt.ylabel("Tanimoto Similarity")
    plt.grid(True, linestyle='--', alpha=0.7)

    plt.ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"\nComparison plot saved as: {output_file}")


def plot_tanimoto_heatmap(similarity_matrix, output_file="tanimoto_heatmap.png"):
    plt.figure(figsize=(10, 8), dpi=200)

    colormap = sns.diverging_palette(300, 145, s=90, l=80, sep=30, center="dark", as_cmap=True)
    ax = sns.heatmap(
        similarity_matrix,
        square=True,
        cmap=colormap,
        vmin=0,
        vmax=1,
        center=0.5,
        xticklabels=True,
        yticklabels=True,
        annot=False,
        cbar_kws={'label': ''}
    )

    ax.invert_yaxis()
    plt.yticks(rotation="horizontal")

    plt.title("Tanimoto Similarity Matrix")
    plt.xlabel("Frame Number")
    plt.ylabel("Frame Number")

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"\nHeatmap saved as: {output_file}")


def analyze_tanimoto_vs_reference(file_path, reference_column=4, step=10):
    try:
        print(f"Reading file: {file_path}")
        df = pd.read_csv(file_path)

        if len(df.columns) <= reference_column:
            raise ValueError(f"El archivo debe tener al menos {reference_column + 1} columnas")

        data = df.iloc[:, reference_column - 1:]

        similarity_series = create_similarity_vector(data, step)
        similarity_matrix = create_similarity_matrix(data, step)

        plot_tanimoto_comparison(similarity_series)
        plot_tanimoto_heatmap(similarity_matrix)
        tanimoto_indices = similarity_series.to_dict()

        print(f"Reference frame: 0")
        print(f"Frames analyzed at {step}-frame intervals")
        print(f"Total frames analyzed: {len(tanimoto_indices)}")

        return tanimoto_indices, similarity_matrix

    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        return None, None


def save_results(tanimoto_indices, output_file="tanimoto_results.txt"):
    try:
        with open(output_file, 'w') as f:
            f.write(f"Tanimoto Analysis Results (10-frame intervals)\n")
            f.write("\n{:<15} {:<10}\n".format("Frame", "Tanimoto"))
            f.write("-" * 25 + "\n")

            for frame in sorted(tanimoto_indices.keys()):
                f.write("{:<15} {:<10}\n".format(frame, round(tanimoto_indices[frame], 4)))

        print(f"\nResults saved to: {output_file}")

    except Exception as e:
        print(f"Error saving results: {str(e)}")


def main():
    reference_column = 4
    step = 10
    tanimoto_indices, similarity_matrix = analyze_tanimoto_vs_reference(FILE_PATH, reference_column, step)

    if tanimoto_indices is not None and similarity_matrix is not None:
        print("\nTanimoto Indices (compared to Frame 0):")
        for frame in sorted(tanimoto_indices.keys()):
            print(f"Frame {frame}: {round(tanimoto_indices[frame], 4)}")
        save_results(tanimoto_indices)


if __name__ == "__main__":
    main()
