from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist

INPUT_CSV = "centered_matrix.csv"
OUTPUT_PDF = "Obese_vs_Diabetic_individual_heatmap_hierarchical.pdf"
COLOR_LIMIT = 2.0
HEATMAP_COLORS = ["#2166AC", "#F7F7F7", "#B2182B"]

CLUSTER_METRIC = "euclidean"
LINKAGE_METHOD = "ward"

OBESE_SAMPLES = [
    "O1", "O2", "O3", "O5", "O6",
    "O7", "O8", "O9", "O10"
]

DIABETIC_SAMPLES = [
    "D1", "D2", "D3", "D5",
    "D6", "D7", "D8S"
]

SAMPLE_ORDER = OBESE_SAMPLES + DIABETIC_SAMPLES

COHORTS = [
    ("Obese", 0, 8),
    ("Diabetic", 9, 15)
]


def orient_linkage_by_score(tree, scores):
    tree = tree.copy()
    n_leaves = len(scores)

    node_scores = {
        i: (float(scores[i]), 1)
        for i in range(n_leaves)
    }

    for row_index in range(tree.shape[0]):
        node_id = n_leaves + row_index

        left_id = int(tree[row_index, 0])
        right_id = int(tree[row_index, 1])

        left_score, left_size = node_scores[left_id]
        right_score, right_size = node_scores[right_id]

        if left_score < right_score:
            tree[row_index, 0], tree[row_index, 1] = (
                tree[row_index, 1],
                tree[row_index, 0]
            )

            left_score, right_score = right_score, left_score
            left_size, right_size = right_size, left_size

        combined_score = (
            left_score * left_size + right_score * right_size
        ) / (left_size + right_size)

        node_scores[node_id] = (
            combined_score,
            left_size + right_size
        )

    return tree


def cluster_order(values, scores=None):
    if values.shape[0] < 2:
        return np.arange(values.shape[0]), None

    distances = pdist(values, metric=CLUSTER_METRIC)

    tree = linkage(
        distances,
        method=LINKAGE_METHOD
    )

    if scores is not None:
        tree = orient_linkage_by_score(tree, scores)

    return leaves_list(tree), tree


def main():
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(INPUT_CSV)

    if not input_path.exists():
        raise FileNotFoundError(f"Cannot find {input_path}")

    df = pd.read_csv(input_path)

    for col in ["gene_symbol", "major_class"]:
        if col not in df.columns:
            raise ValueError(f"Missing required annotation column: {col}")

    missing_samples = [
        sample for sample in SAMPLE_ORDER
        if sample not in df.columns
    ]

    if missing_samples:
        raise ValueError(
            "Missing sample columns: " + ", ".join(missing_samples)
        )

    matrix = df[SAMPLE_ORDER].apply(
        pd.to_numeric,
        errors="coerce"
    ).to_numpy(dtype=float)

    if not np.isfinite(matrix).all():
        raise ValueError(
            "All Obese and Diabetic sample columns must contain finite numeric values."
        )

    row_labels = [
        f"{symbol}  [{str(category).replace(' channels', '').replace('; ', ' + ')}]"
        for symbol, category in zip(
            df["gene_symbol"],
            df["major_class"]
        )
    ]

    row_scores = matrix.mean(axis=1)

    row_order, row_tree = cluster_order(
        matrix,
        scores=row_scores
    )

    clustered_row_matrix = matrix[row_order, :]

    column_order = []
    cohort_trees = []

    for cohort_name, start, end in COHORTS:
        cohort_indices = np.arange(start, end + 1)

        cohort_matrix = clustered_row_matrix[:, cohort_indices].T

        local_order, cohort_tree = cluster_order(cohort_matrix)

        ordered_indices = cohort_indices[local_order]

        column_order.extend(ordered_indices.tolist())

        cohort_trees.append(
            (cohort_name, cohort_tree, len(ordered_indices))
        )

    column_order = np.array(column_order)

    matrix = matrix[np.ix_(row_order, column_order)]

    row_labels = np.array(row_labels)[row_order].tolist()
    sample_order = np.array(SAMPLE_ORDER)[column_order].tolist()

    cohort_codes_original = np.array(
        [0] * len(OBESE_SAMPLES) +
        [1] * len(DIABETIC_SAMPLES)
    )

    cohort_codes = cohort_codes_original[column_order]

    n_genes = len(row_labels)
    n_samples = len(sample_order)

    cohort_sizes = [
        len(OBESE_SAMPLES),
        len(DIABETIC_SAMPLES)
    ]

    boundaries = np.cumsum(cohort_sizes)[:-1]

    fig_height = max(20, 0.22 * n_genes + 6)

    fig = plt.figure(figsize=(18, fig_height))

    gs = fig.add_gridspec(
        4,
        4,
        height_ratios=[0.55, 0.9, 0.18, 14],
        width_ratios=[0.18, 0.28, 1.0, 0.035],
        hspace=0.06,
        wspace=0.025
    )

    fig.subplots_adjust(
        left=0.10,
        right=0.96,
        top=0.975,
        bottom=0.05
    )

    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")

    title_ax.text(
        0.5,
        0.70,
        "Person-by-Person Ion-Channel / Transporter Expression Heatmap: Obese vs Diabetic",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold"
    )

    title_ax.text(
        0.5,
        0.20,
        "Genes hierarchically clustered by Obese + Diabetic expression pattern and oriented from higher to lower mean centered expression; samples clustered within cohort",
        ha="center",
        va="center",
        fontsize=9
    )

    dendrogram_grid = gs[1, 2].subgridspec(
        1,
        2,
        width_ratios=cohort_sizes,
        wspace=0.03
    )

    for i, (cohort_name, cohort_tree, cohort_size) in enumerate(cohort_trees):
        dend_ax = fig.add_subplot(dendrogram_grid[0, i])

        if cohort_tree is not None:
            dendrogram(
                cohort_tree,
                orientation="top",
                no_labels=True,
                color_threshold=0,
                above_threshold_color="black",
                ax=dend_ax
            )

        dend_ax.set_xticks([])
        dend_ax.set_yticks([])

        for spine in dend_ax.spines.values():
            spine.set_visible(False)

    row_dend_ax = fig.add_subplot(gs[3, 0])

    if row_tree is not None:
        dendrogram(
            row_tree,
            orientation="left",
            no_labels=True,
            color_threshold=0,
            above_threshold_color="black",
            ax=row_dend_ax
        )

        row_dend_ax.set_ylim(10 * n_genes, 0)

    row_dend_ax.set_xticks([])
    row_dend_ax.set_yticks([])

    for spine in row_dend_ax.spines.values():
        spine.set_visible(False)

    labels_ax = fig.add_subplot(gs[3, 1])
    labels_ax.set_xlim(0, 1)
    labels_ax.set_ylim(n_genes - 0.5, -0.5)
    labels_ax.axis("off")

    font_size = max(5.2, min(8.0, 190 / max(n_genes, 1)))

    for row_index, label in enumerate(row_labels):
        labels_ax.text(
            0.99,
            row_index,
            label,
            ha="right",
            va="center",
            fontsize=font_size
        )

    strip_ax = fig.add_subplot(gs[2, 2])

    strip_ax.imshow(
        cohort_codes.reshape(1, -1),
        aspect="auto",
        cmap=plt.get_cmap("tab10", 2),
        vmin=0,
        vmax=1
    )

    strip_ax.set_xlim(-0.5, n_samples - 0.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])

    for boundary in boundaries:
        strip_ax.axvline(
            boundary - 0.5,
            color="black",
            linewidth=1.0
        )

    start = 0

    for cohort_name, cohort_size in zip(
        ["Obese", "Diabetic"],
        cohort_sizes
    ):
        middle = start + (cohort_size - 1) / 2

        strip_ax.text(
            middle,
            1.55,
            cohort_name,
            transform=strip_ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold"
        )

        start += cohort_size

    for spine in strip_ax.spines.values():
        spine.set_visible(False)

    heat_ax = fig.add_subplot(gs[3, 2])

    cmap = LinearSegmentedColormap.from_list(
        "blue_white_red",
        HEATMAP_COLORS,
        N=256
    )

    norm = TwoSlopeNorm(
        vmin=-COLOR_LIMIT,
        vcenter=0,
        vmax=COLOR_LIMIT
    )

    image = heat_ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        norm=norm
    )

    heat_ax.set_xticks(np.arange(n_samples))

    heat_ax.set_xticklabels(
        sample_order,
        rotation=90,
        fontsize=8
    )

    heat_ax.tick_params(
        axis="x",
        length=0,
        pad=3
    )

    heat_ax.set_yticks([])
    heat_ax.set_xlim(-0.5, n_samples - 0.5)
    heat_ax.set_ylim(n_genes - 0.5, -0.5)

    for boundary in boundaries:
        heat_ax.axvline(
            boundary - 0.5,
            color="black",
            linewidth=1.0
        )

    cbar_ax = fig.add_subplot(gs[3, 3])

    cbar = fig.colorbar(
        image,
        cax=cbar_ax
    )

    cbar.set_label(
        "Centered expression\n(log2-CPM deviation from the Normal mean)",
        fontsize=9
    )

    cbar.ax.tick_params(labelsize=8)

    fig.text(
        0.61,
        0.012,
        "Blue = lower than the Normal-group mean; white = near the Normal-group mean; "
        "red = higher. "
        f"Color range: −{COLOR_LIMIT:g} to +{COLOR_LIMIT:g}; values beyond this range are color-saturated.",
        ha="center",
        va="bottom",
        fontsize=8
    )

    with PdfPages(OUTPUT_PDF) as pdf:
        pdf.savefig(
            fig,
            bbox_inches="tight",
            pad_inches=0.15
        )

    plt.close(fig)


if __name__ == "__main__":
    main()