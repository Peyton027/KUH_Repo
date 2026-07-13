library(ggplot2)
library(dplyr)
library(ggrepel)

rna <- read.delim("AP_NvHF.tabular")
protein <- read.csv("HfvN_prot.csv", check.names = TRUE)

protein <- protein[, !grepl("^X$|^Unnamed", names(protein))]

rna$Gene <- toupper(trimws(rna$Gene))
protein$Gene <- toupper(trimws(protein$Gene))

rna$RNA_log2FC <- as.numeric(rna$RNA_log2FC)
protein$Protein_log2FC <- as.numeric(protein$Log_FC)

rna <- rna %>%
  filter(!is.na(Gene), !is.na(RNA_log2FC)) %>%
  group_by(Gene) %>%
  summarise(RNA_log2FC = median(RNA_log2FC, na.rm = TRUE))

protein <- protein %>%
  filter(!is.na(Gene), !is.na(Protein_log2FC)) %>%
  group_by(Gene) %>%
  summarise(Protein_log2FC = median(Protein_log2FC, na.rm = TRUE))

combined <- inner_join(rna, protein, by = "Gene")

combined <- combined %>%
  filter(
    is.finite(RNA_log2FC),
    is.finite(Protein_log2FC)
  )

combined$Behavior <- ifelse(
  combined$RNA_log2FC > 0 & combined$Protein_log2FC > 0,
  "Up in both",
  ifelse(
    combined$RNA_log2FC < 0 & combined$Protein_log2FC < 0,
    "Down in both",
    "Opposite directions"
  )
)

spearman_test <- cor.test(
  combined$RNA_log2FC,
  combined$Protein_log2FC,
  method = "spearman",
  exact = FALSE
)

rho <- round(spearman_test$estimate, 2)
p_value <- spearman_test$p.value

if (p_value < 0.001) {
  p_text <- "p < 0.001"
} else {
  p_text <- paste0("p = ", round(p_value, 3))
}

cor_text <- paste0(
  "Spearman rho = ", rho,
  "\n", p_text,
  "\nn = ", nrow(combined)
)

genes_to_label <- c(
  "ATP1A1",
  "ATP1A2",
  "ATP1A3",
  "ATP1A4",
  "ATP1B1",
  "ATP1B2",
  "ATP1B3",
  "ATP1B4",
  "SCN5A"
)

label_data <- combined %>%
  filter(Gene %in% genes_to_label)

plot_limit <- max(
  abs(combined$RNA_log2FC),
  abs(combined$Protein_log2FC)
)

plot_limit <- plot_limit * 1.15

p <- ggplot(
  combined,
  aes(
    x = RNA_log2FC,
    y = Protein_log2FC,
    color = Behavior
  )
) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    data = label_data,
    aes(label = Gene),
    color = "black",
    size = 4,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = -plot_limit * 0.95,
    y = plot_limit * 0.95,
    label = cor_text,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "Up in both" = "red3",
      "Down in both" = "blue3",
      "Opposite directions" = "orange"
    )
  ) +
  xlim(-plot_limit, plot_limit) +
  ylim(-plot_limit, plot_limit) +
  coord_equal() +
  labs(
    title = "Transcriptomic-Proteomic Concordance",
    subtitle = "Heart failure vs normal",
    x = "RNA log2 fold change",
    y = "Protein log2 fold change",
    color = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "bottom"
  )

print(p)

ggsave(
  "RNA_Protein_Concordance.pdf",
  p,
  width = 7,
  height = 7
)

write.csv(
  combined,
  "RNA_Protein_Concordance_Data.csv",
  row.names = FALSE
)
