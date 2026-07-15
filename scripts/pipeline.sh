#A WIP rudimentary shell script for workflow automation 


#!/bin/bash

set -e

PROJECT_DIR="$(pwd)"
DATA_DIR="$PROJECT_DIR/data"
SCRIPT_DIR="$PROJECT_DIR/scripts"
RESULTS_DIR="$PROJECT_DIR/results"
FIGURE_DIR="$RESULTS_DIR/figures"

RNA_COUNTS="$DATA_DIR/rna_counts.csv"
RNA_METADATA="$DATA_DIR/sample_metadata.tsv"
PROTEOMICS="$DATA_DIR/proteomics.csv"
GENE_LIST="$DATA_DIR/ion_handling_genes.csv"

mkdir -p "$RESULTS_DIR"
mkdir -p "$FIGURE_DIR"

echo "Starting cardiac ion-handling analysis"

echo "Checking project files"
ls "$DATA_DIR"
ls "$SCRIPT_DIR"

echo "Preparing transcriptomic analysis"
Rscript "$SCRIPT_DIR/pca_analysis.R" "$RNA_COUNTS" "$RNA_METADATA"
Rscript "$SCRIPT_DIR/heatmap_analysis.R" "$RNA_COUNTS" "$GENE_LIST"

echo "Preparing differential-expression analysis"
Rscript "$SCRIPT_DIR/deseq2_analysis.R" "$RNA_COUNTS" "$RNA_METADATA"

echo "Preparing proteomics analysis"
python "$SCRIPT_DIR/proteomics_analysis.py" "$PROTEOMICS" "$GENE_LIST"

echo "Additional analysis steps remain in development"
