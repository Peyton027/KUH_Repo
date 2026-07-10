{r}

library(factoextra)
library(MASS)
library(dplyr)
library(ggplot2)

# load all genes
all_rna <- read.csv("RNA seq raw file.csv")

# load ion-channel genes
ion <- read.csv("RNA_ion_channels.csv")

# save gene id as row names
rownames(all_rna) <- all_rna$gene_id
rownames(ion) <- ion$gene_id

# get Normal and HF sample names from ion-channel file
sample_names <- colnames(ion)[-1]

# keep only the same Normal and HF sample columns
all_rna <- all_rna[, sample_names]
ion <- ion[, sample_names]

# find library size using all genes
smpl_size <- colSums(all_rna)

# normalize ion-channel genes
rna_cpm <- t(t(ion) / smpl_size * 1000000)

# log2 transform
rna_log <- log2(rna_cpm + 1)

# get rid of 0 variance genes
genes <- apply(rna_log, 1, var, na.rm = TRUE) != 0
rna_log <- rna_log[genes, ]

dim(rna_log)
head(rna_log)

# transpose matrix
pca <- prcomp(
  t(rna_log),
  scale. = TRUE
)

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


