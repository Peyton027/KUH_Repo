import pandas as pd


raw = pd.read_csv("RNA-seq_raw.csv")
normal_samples = ["N1", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"]

hf_samples = [
    "HF1", "HF2", "HF3", "HF4", "HF5",
    "HF6", "HF7", "HF8", "HF9", "HF10"
]

samples = normal_samples + hf_samples


# hardcoded target genes

target_genes = [
    # Phase 0: Sodium current
    "SCN10A", "SCN1B", "SCN5A", "SCN2B", "SCN3B",
    "SCN4B", "SCN11A", "SCN8A",

    # Phase 1: Early repolarization
    "KCNA5", "KCNC1", "KCNAB1", "KCNAB2",
    "KCND2", "KCND3",
    "KCNJ11",
    "KCNA4", "KCNA7", "KCNC4",
    "DPP10", "DPP6", "KCNIP2",

    # Phase 2: Plateau / calcium entry
    "CACNA1B", "CACNA1C",
    "CACNA2D1", "CACNB2",
    "PKD2", "TRPC1", "TRPM4",

    # All-phase calcium handling
    "SLC8A1", "SLC8A2", "SLC8A3",
    "ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4",
    "ATP2A1", "ATP2A2", "ATP2A3",
    "SLC8B1",

    # Phase 3: Delayed rectifier / inward rectifier K+ currents
    "KCNH2",
    "KCNE1", "KCNE2", "KCNE3", "KCNE4", "KCNE5",
    "KCNQ1",
    "KCNJ12", "KCNJ14", "KCNJ2", "KCNJ4",

    # Phase 4: Late repolarization / pacemaker / K+ currents
    "CACNA1G", "CACNA1H",
    "HCN1", "HCN2", "HCN3", "HCN4",
    "KCNJ3", "KCNJ5",
    "ABCC8", "KCNJ8",
    "KCNK1", "KCNK3", "KCNK6",
    "KCNK17", "KCNK9",
    "KCNMA1", "KCNMB1",
    "KCNN1", "KCNN3",

    # All-phase Na+/K+-ATPase accessory subunits
    "ATP1B1", "ATP1B2", "ATP1B3", "ATP1B4"
]


required_columns = [
    "gene_id",
    "gene_name",
    "gene_description",
    "gene_biotype"
] + samples

missing_columns = [
    column for column in required_columns
    if column not in raw.columns
]

if missing_columns:
    raise ValueError(
        "RNA-seq_raw.csv is missing these required columns:\n"
        + ", ".join(missing_columns)
    )


raw["gene_name"] = (
    raw["gene_name"]
    .astype(str)
    .str.strip()
    .str.upper()
)


selected = raw[
    raw["gene_name"].isin(target_genes)
].copy()

selected["gene_name"] = pd.Categorical(
    selected["gene_name"],
    categories=target_genes,
    ordered=True
)

selected = selected.sort_values("gene_name")

# QC - for missing genes

found_genes = set(selected["gene_name"].astype(str))

missing_target_genes = [
    gene for gene in target_genes
    if gene not in found_genes
]

if missing_target_genes:
    print("\nWARNING: These target genes were not found in RNA-seq_raw.csv:")
    print(", ".join(missing_target_genes))
else:
    print("\nAll target genes were found in RNA-seq_raw.csv.")

print(f"\nTarget genes requested: {len(target_genes)}")
print(f"Target genes found: {len(found_genes)}")


# export cnts for deseq2

counts_matrix = selected[
    ["gene_id"] + samples
]

counts_matrix.to_csv(
    "Galaxy_DESeq2_Normal_vs_HF_counts_matrix.tsv",
    sep="\t",
    index=False
)

# save metadata

metadata = pd.DataFrame({
    "sample_id": samples,
    "condition": (
        ["Normal"] * len(normal_samples)
        + ["HF"] * len(hf_samples)
    )
})

metadata.to_csv(
    "Galaxy_DESeq2_Normal_vs_HF_sample_metadata.tsv",
    sep="\t",
    index=False
)

# keeps gene annotation

annotation = selected[
    [
        "gene_id",
        "gene_name",
        "gene_description",
        "gene_biotype"
    ]
].copy()

annotation["gene_name"] = annotation["gene_name"].astype(str)

annotation.to_csv(
    "Normal_vs_HF_target_gene_annotation.tsv",
    sep="\t",
    index=False
)

# EXPORT GENE-AVAILABILITY AUDIT

audit = pd.DataFrame({
    "gene_name": target_genes,
    "status_in_raw_matrix": [
        "Found" if gene in found_genes else "Missing"
        for gene in target_genes
    ]
})

audit.to_csv(
    "Normal_vs_HF_target_gene_audit.tsv",
    sep="\t",
    index=False
)

print("\nFiles created:")
print("1. Galaxy_DESeq2_Normal_vs_HF_counts_matrix.tsv")
print("2. Galaxy_DESeq2_Normal_vs_HF_sample_metadata.tsv")
print("3. Normal_vs_HF_target_gene_annotation.tsv")
print("4. Normal_vs_HF_target_gene_audit.tsv")