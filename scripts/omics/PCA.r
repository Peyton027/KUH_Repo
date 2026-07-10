library(factoextra)
library(MASS)
library(dplyr)
library(ggplot2)

# load data 
raw_rna <- read.csv("RNA_ion_channels.csv")
rownames(raw_rna) <- raw_rna$gene_id
raw_rna <- raw_rna[, -1]

head(raw_rna)
dim(raw_rna)


#keep only genes with nonzero variance
genes <-apply(raw_rna,1,var, na.rm=T)!=0
rna <-rna[genes,]
dim(rna)
head(rna)



#normalize

smpl_size <-colSums(rna)
rna_cpm <-t(t(rna)/smpl_size *1000000)

#log2 transform 
rna_log <- log2(rna_cpm + 1)

#transpose matrix 
pca <- prcomp(t(rna_log), scale. = TRUE)

#add condition class 
sample_names <- rownames(pca$x)

condition <- ifelse(
  grepl("^N", sample_names),
  "Normal",
  "HF"
)

#scree plot 
fviz_eig(pca,addlabels = TRUE)

#pca plot
fviz_pca_ind(
  pca,
  geom.ind = "point",
  habillage = condition,
  addEllipses = TRUE,
  repel = T
)
