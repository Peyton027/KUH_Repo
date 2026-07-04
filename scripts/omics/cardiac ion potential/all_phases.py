from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import pdist

DESEQ_FILE = "AP_NvHF.tabular"
ANNOTATION_FILE = "ap_channel_manuscript_table_all_phases_updated.csv"

PHASE_TO_PLOT = "All phases"
OUTPUT_PDF = "Normal_vs_HF_AllPhases_mirrored_DESeq2_heatmap.pdf"

COLOR_LIMIT = 2.0
HEATMAP_COLORS = ["#2166AC", "#F7F7F7", "#B2182B"]

GENE_TO_ENSEMBL = {
    "SLC8A1": "ENSG00000183023",
    "SLC8A2": "ENSG00000118160",
    "SLC8A3": "ENSG00000100678",
    "SLC8B1": "ENSG00000089060",
    "ATP2B1": "ENSG00000070961",
    "ATP2B2": "ENSG00000157087",
    "ATP2B3": "ENSG00000067842",
    "ATP2B4": "ENSG00000058668",
    "ATP2A1": "ENSG00000196296",
    "ATP2A2": "ENSG00000174437",
    "ATP2A3": "ENSG00000074370",
    "KCNK1": "ENSG00000135750",
    "KCNK3": "ENSG00000171303",
    "KCNK6": "ENSG00000099337",
    "KCNK17": "ENSG00000124780",
    "KCNK9": "ENSG00000169427",
    "ATP1B1": "ENSG00000143153",
    "ATP1B2": "ENSG00000129244",
    "ATP1B3": "ENSG00000069849",
    "ATP1B4": "ENSG00000101892"
}


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


def cluster_rows(matrix, orientation_scores):
    if matrix.shape[0] < 2:
        return np.arange(matrix.shape[0]), None

    tree = linkage(
        pdist(matrix, metric="euclidean"),
        method="ward"
    )

    tree = orient_linkage_by_score(tree, orientation_scores)

    return leaves_list(tree), tree


def get_fdr_stars(padj):
    if pd.isna(padj):
        return ""
    elif padj < 0.0001:
        return " ****"
    elif padj < 0.001:
        return " ***"
    elif padj < 0.01:
        return " **"
    elif padj < 0.05:
        return " *"
    return ""


def main():
    for file_path in [DESEQ_FILE, ANNOTATION_FILE]:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Cannot find: {file_path}")

    annotation_df = pd.read_csv(ANNOTATION_FILE)

    for col in ["phase_group", "gene", "current"]:
        if col not in annotation_df.columns:
            raise ValueError(f"Missing required annotation column: {col}")

    phase_df = annotation_df.loc[
        annotation_df["phase_group"].astype(str).str.strip().eq(PHASE_TO_PLOT),
        ["gene", "current"]
    ].copy()

    if phase_df.empty:
        raise ValueError(f"No genes found for {PHASE_TO_PLOT}")

    phase_df["gene"] = (
        phase_df["gene"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    phase_df = phase_df.drop_duplicates(subset="gene", keep="first")
    phase_df["gene_id"] = phase_df["gene"].map(GENE_TO_ENSEMBL)

    unmapped_genes = phase_df.loc[
        phase_df["gene_id"].isna(),
        "gene"
    ].tolist()

    if unmapped_genes:
        print("Warning: No Ensembl mapping for: " + ", ".join(unmapped_genes))

    phase_df = phase_df.dropna(subset=["gene_id"]).copy()

    deseq_columns = [
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    ]

    deseq_df = pd.read_csv(
        DESEQ_FILE,
        sep="\t",
        header=None,
        names=deseq_columns
    )

    deseq_df["gene_id"] = (
        deseq_df["gene_id"]
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
    )

    deseq_df = deseq_df.drop_duplicates(subset="gene_id", keep="first")

    df = phase_df.merge(
        deseq_df[["gene_id", "log2FoldChange", "pvalue", "padj"]],
        on="gene_id",
        how="left"
    )

    missing_deseq = df.loc[df["log2FoldChange"].isna(), "gene"].tolist()

    if missing_deseq:
        print("Warning: No DESeq2 value found for: " + ", ".join(missing_deseq))

    df = df.dropna(subset=["log2FoldChange"]).copy()

    if df.empty:
        raise ValueError("No usable DESeq2 values found for the All phases group.")

    hf_values = df["log2FoldChange"].to_numpy(dtype=float)

    matrix = np.column_stack([
        -hf_values,
        hf_values
    ])

    row_labels = [
        f"{gene}  [{current}]{get_fdr_stars(padj)}"
        for gene, current, padj in zip(
            df["gene"],
            df["current"],
            df["padj"]
        )
    ]

    row_order, row_tree = cluster_rows(matrix, hf_values)

    matrix = matrix[row_order, :]
    row_labels = np.array(row_labels)[row_order].tolist()

    n_genes = len(row_labels)
    group_labels = ["Normal", "Heart Failure"]

    fig_height = max(8, 0.55 * n_genes + 5)
    fig = plt.figure(figsize=(12, fig_height))

    gs = fig.add_gridspec(
        4,
        4,
        height_ratios=[0.60, 0.12, 0.18, 12],
        width_ratios=[0.20, 0.42, 1.0, 0.06],
        hspace=0.06,
        wspace=0.03
    )

    fig.subplots_adjust(
        left=0.10,
        right=0.94,
        top=0.97,
        bottom=0.08
    )

    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")

    title_ax.text(
        0.5,
        0.70,
        "All-Phases Group Differential-Expression Heatmap: Normal vs Heart Failure",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold"
    )

    title_ax.text(
        0.5,
        0.20,
        "Only genes labeled All phases; Heart Failure = DESeq2 log2FC, Normal = inverse",
        ha="center",
        va="center",
        fontsize=9
    )

    strip_ax = fig.add_subplot(gs[2, 2])

    strip_ax.imshow(
        np.array([[0, 1]]),
        aspect="auto",
        cmap=plt.get_cmap("tab10", 2),
        vmin=0,
        vmax=1
    )

    strip_ax.set_xlim(-0.5, 1.5)
    strip_ax.set_xticks([])
    strip_ax.set_yticks([])

    for index, label in enumerate(group_labels):
        strip_ax.text(
            index,
            1.55,
            label,
            transform=strip_ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold"
        )

    strip_ax.axvline(0.5, color="black", linewidth=1.0)

    for spine in strip_ax.spines.values():
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

    font_size = max(7.0, min(10.0, 150 / max(n_genes, 1)))

    for row_index, label in enumerate(row_labels):
        labels_ax.text(
            0.99,
            row_index,
            label,
            ha="right",
            va="center",
            fontsize=font_size
        )

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

    heat_ax.set_xticks([0, 1])
    heat_ax.set_xticklabels(group_labels, fontsize=10)
    heat_ax.tick_params(axis="x", length=0, pad=8)
    heat_ax.set_yticks([])
    heat_ax.set_xlim(-0.5, 1.5)
    heat_ax.set_ylim(n_genes - 0.5, -0.5)
    heat_ax.axvline(0.5, color="black", linewidth=1.0)

    cbar_ax = fig.add_subplot(gs[3, 3])
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label("Mirrored DESeq2 log2FC", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.text(
        0.60,
        0.018,
        "Heart Failure = DESeq2 log2FC. Normal = inverse of that value. "
        "This is a mirrored differential-expression display, not raw expression. "
        "* padj < 0.05; ** padj < 0.01; *** padj < 0.001; **** padj < 0.0001.",
        ha="center",
        va="bottom",
        fontsize=8
    )

    with PdfPages(OUTPUT_PDF) as pdf:
        pdf.savefig(fig, bbox_inches="tight", pad_inches=0.15)

    plt.close(fig)

    print(
        f"Saved {PHASE_TO_PLOT} heatmap with "
        f"{n_genes} genes to: {OUTPUT_PDF}"
    )


if __name__ == "__main__":
    main()