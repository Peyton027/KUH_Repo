library(factoextra)
library(MASS)
library(dplyr)
library(ggplot2)
library(edgeR)

# load all genes
all_rna <- read.csv(
  "RNA seq raw file.csv",
  check.names = FALSE
)

# load ion-channel genes
ion <- read.csv(
  "RNA_ion_channels.csv",
  check.names = FALSE
)

# save gene id as row names
rownames(all_rna) <- all_rna$gene_id
rownames(ion) <- ion$gene_id

# get Normal and HF sample names from ion-channel file
sample_names <- colnames(ion)[-1]

# keep only the same Normal and HF sample columns
all_rna <- all_rna[, sample_names]
ion <- ion[, sample_names]

# make sure counts are numeric
all_rna[] <- lapply(all_rna, as.numeric)
ion[] <- lapply(ion, as.numeric)

# create edgeR object using all genes
dge <- DGEList(counts = all_rna)

# calculate TMM normalization factors using all genes
dge <- calcNormFactors(
  dge,
  method = "TMM"
)

# calculate TMM-normalized log2 CPM for all genes
rna_log_all <- cpm(
  dge,
  log = TRUE,
  prior.count = 2
)

# find ion-channel genes present in the complete RNA-seq file
ion_genes <- intersect(
  rownames(ion),
  rownames(rna_log_all)
)

# keep only ion-channel genes
rna_log <- rna_log_all[ion_genes, ]

# get rid of 0 variance genes
genes <- apply(
  rna_log,
  1,
  var,
  na.rm = TRUE
) != 0

rna_log <- rna_log[genes, ]

dim(rna_log)
head(rna_log)

# transpose matrix and perform PCA
pca <- prcomp(t(rna_log),scale. = TRUE)

# add condition class
sample_names <- rownames(pca$x)

condition <- ifelse(
  grepl("^N", sample_names),
  "Normal",
  "HF"
)

# scree plot
fviz_eig(
  pca,
  addlabels = TRUE
)

# PCA plot
fviz_pca_ind(
  pca,
  geom.ind = "point",
  habillage = condition,
  addEllipses = TRUE,
  repel = TRUE
)
