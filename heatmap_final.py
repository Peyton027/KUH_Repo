import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

INPUT_FILE = "ion_channel_heatmap_log2FC_annotated.csv"
OUTPUT_FILE = "Ranked_6column_Heatmap_compact.pdf"

COLOR_MIN = -2
COLOR_MAX = 2

comparisons = [
    {
        "condition": "Obese",
        "lfc": "Obese_vs_Normal_log2FoldChange"
    },
    {
        "condition": "Diabetes",
        "lfc": "Diabetic_vs_Normal_log2FoldChange"
    },
    {
        "condition": "HF",
        "lfc": "HF_vs_Normal_log2FoldChange"
    }
]

df = pd.read_csv(INPUT_FILE)

sorted_data = []

for comparison in comparisons:
    sub = df[[comparison["lfc"]]].dropna().copy()
    sub = sub.sort_values(comparison["lfc"], ascending=False)
    sorted_data.append(sub.reset_index(drop=True))

max_rows = max(len(sub) for sub in sorted_data)

heatmap_values = np.full((max_rows, 6), np.nan)

for i, sub in enumerate(sorted_data):
    lfc_values = sub[comparisons[i]["lfc"]].values

    heatmap_values[:len(lfc_values), i * 2] = -lfc_values
    heatmap_values[:len(lfc_values), i * 2 + 1] = lfc_values

cmap = LinearSegmentedColormap.from_list(
    "stronger_blue_white_red",
    [
        (0.00, "#001A66"),
        (0.18, "#0047AB"),
        (0.35, "#4F9BD6"),
        (0.47, "#D8ECF8"),
        (0.50, "#FFFFFF"),
        (0.53, "#FDE0DD"),
        (0.65, "#F28E8B"),
        (0.82, "#E53935"),
        (1.00, "#8B0000")
    ]
)

cmap.set_bad("white")

fig_height = max(8, max_rows * 0.055)

fig, ax = plt.subplots(
    figsize=(5.6, fig_height),
    facecolor="white"
)

im = ax.imshow(
    heatmap_values,
    aspect="auto",
    cmap=cmap,
    vmin=COLOR_MIN,
    vmax=COLOR_MAX,
    interpolation="none"
)

ax.set_xticks([0, 1, 2, 3, 4, 5])

ax.set_xticklabels(
    [
        "Normal-Obese",
        "Obese",
        "Normal-Db",
        "Diabetes",
        "Normal-HF",
        "HF"
    ],
    fontsize=9,
    rotation=35,
    ha="right"
)

ax.set_yticks([])

ax.tick_params(
    axis="x",
    length=0,
    pad=6
)

ax.vlines(
    [1.5, 3.5],
    -0.5,
    max_rows - 0.5,
    color="black",
    linewidth=0.8
)

for spine in ax.spines.values():
    spine.set_visible(False)

cbar = plt.colorbar(
    im,
    ax=ax,
    pad=0.025,
    fraction=0.035,
    ticks=[-2, -1, 0, 1, 2]
)

cbar.set_label(
    "log2FC",
    fontsize=9
)

cbar.ax.tick_params(
    labelsize=8,
    length=0
)

fig.suptitle(
    "Ion Channel Gene Heatmap: Ranked Differential Expression",
    fontsize=12,
    y=0.985
)

fig.subplots_adjust(
    left=0.04,
    right=0.88,
    top=0.95,
    bottom=0.08
)

plt.savefig(
    OUTPUT_FILE,
    bbox_inches="tight"
)

plt.show()