# Integrative Dual-Omics Analysis Reveals Ion Channel Remodeling Across the Cardiomyocyte Action Potential

**Author:** Peyton Chiu  
**Affiliation:** The University of Texas Health Science Center at San Antonio  
**Home Institution:** University of California San Diego, B.S. Bioinformatics  
**Presentation:** NIDDK KUH Undergraduate Research Symposium, University of Pittsburgh, July 23–25, 2026

## Project Overview

Heart failure is associated with structural and electrical remodeling that can disrupt the cardiac action potential and increase susceptibility to arrhythmias. This project integrates transcriptomic and proteomic data from human left-ventricular tissue to examine how cardiac ion channels, pumps, transporters, and associated regulatory genes are altered in heart failure with reduced ejection fraction (HFrEF).

## Study Design

Human left-ventricular samples from non-failing controls and patients with HFrEF were analyzed using:

- RNA sequencing
- LC-MS/MS proteomics
- Targeted analysis of cardiac ion-handling genes
- Literature-based assignment of genes to cardiac action-potential phases

The RNA-seq and proteomics datasets were processed separately and later integrated through RNA–protein concordance analysis.

## Analysis Workflow

<p align="center">
  <img src="figures/KUH%20workflow.png"
       alt="Dual-omics analysis workflow for cardiac ion-handling genes"
       width="1000">
</p>

## Cardiac Action-Potential Groups

Genes were assigned according to their primary physiological roles:

- **Phase 0:** Rapid sodium-dependent depolarization
- **Phase 1:** Initial repolarization
- **Phase 2:** Calcium-dependent plateau phase
- **Phase 3:** Potassium-dependent repolarization
- **Phase 4:** Resting membrane potential and pacemaker activity
- **All phases:** Ion homeostasis, calcium cycling, pumps, exchangers, and regulatory proteins

## Main Analyses

### Transcriptomics

- Cardiac ion-handling genes were extracted from the RNA-seq count matrix.
- Expression values were normalized for targeted PCA and heatmap visualization.
- DESeq2 was used to estimate differential expression between normal-control and HFrEF samples.
- Genes were visualized by cardiac action-potential phase.

### Proteomics

- Protein abundance data were filtered to retain matched cardiac ion-handling proteins.
- Protein-level fold changes were calculated for HFrEF relative to normal controls.
- Detected proteins were mapped to the same cardiac action-potential framework used for the transcriptomic analysis.

### Multi-Omic Integration

Matched RNA and protein log2 fold changes were compared to identify:

- Concordant upregulation
- Concordant downregulation
- Discordant RNA and protein changes
- Genes detected only at the transcript level

## Key Findings

- Targeted PCA showed modest separation between normal-control and HFrEF samples.
- The cohort-centered heatmap showed distinct cardiac ion-handling expression patterns between the two groups.
- **ATP1A1, ATP1A2, and ATP1A3**, which encode Na+/K+-ATPase alpha subunits, were consistently downregulated at both the transcript and protein levels.
- **SCN5A**, which encodes the primary cardiac sodium channel involved in Phase 0 depolarization, was also downregulated.
- **ATP2A2/SERCA2** and **PLN** were reduced, indicating remodeling of sarcoplasmic-reticulum calcium cycling.
- Increased expression of **TRPC1, PKD2, TRPM4**, and calcium-activated potassium-channel genes suggests additional remodeling of calcium-responsive membrane pathways.
- Overall, heart failure was associated with coordinated changes in sodium handling, calcium cycling, and cardiac excitability.

## Repository Structure

```text
.
├── README.md
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
├── scripts/
│   ├── transcriptomics/
│   ├── proteomics/
│   ├── pca/
│   ├── heatmaps/
│   ├── action_potential_mapping/
│   └── concordance/
├── figures/
│   ├── poster_figures/
│   └── supplementary_figures/
├── tables/
│   ├── gene_annotations/
│   └── differential_expression/
└── poster/
    ├── abstract/
    └── final_poster/
```

The exact folder names may be adjusted to match the final repository organization.

## Software

The analysis uses R and Python. Major tools and packages include:

### R

- DESeq2 (Galaxy)
- edgeR
- ggplot2
- factoextra
- dplyr
- ggrepel

### Python

- pandas
- NumPy
- Matplotlib
- SciPy
- openpyxl

Package requirements may vary between scripts.

## Running the Analysis

1. Clone the repository:

```bash
git clone <REPOSITORY_URL>
cd <REPOSITORY_NAME>
```

2. Place the required input files in the appropriate `data/` folders.

3. Confirm that file paths and sample names match those specified in each script.

4. Run the transcriptomic, proteomic, visualization, and concordance scripts from their respective folders.

Example:

```bash
python scripts/proteomics/proteomics_heatmap.py
```

or in R:

```r
source("scripts/transcriptomics/pca_analysis.R")
```

5. Generated figures should be saved to the `figures/` directory.


## Interpretation and Limitations

This study evaluates gene and protein abundance rather than direct ion-channel activity. The results support altered sodium and calcium handling in HFrEF, but intracellular ion concentrations, membrane currents, and action-potential properties were not directly measured. Functional studies are needed to determine whether the observed expression changes produce the predicted electrophysiological effects.

## Acknowledgments

I thank Dr. Muniswamy Madesh and members of the Madesh Laboratory for their mentorship and support. I also thank the UT Health San Antonio SPURRS program and the NIDDK KUH Summer Research Program for providing this research opportunity.

## Funding

This work was supported by the National Institute of Diabetes and Digestive and Kidney Diseases Summer Research Experience Program under grant **R25DK119180**.

## Contact

**Peyton Chiu**  
B.S. Bioinformatics, University of California San Diego  
Expected Graduation: March 2027  
Email: plchiu@ucsd.edu

## Citation

When referencing this repository, please cite the poster:

> Chiu, P. *Integrative Dual-Omics Analysis Reveals Ion Channel Remodeling Across the Cardiomyocyte Action Potential.* NIDDK KUH Undergraduate Research Symposium, University of Pittsburgh, 2026.
