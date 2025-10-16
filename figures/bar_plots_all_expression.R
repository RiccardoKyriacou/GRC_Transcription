library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)

box_theme <- theme_bw() + 
  theme(plot.title = element_text(color="black", size=15),
        panel.background = element_rect(fill="white"),
        axis.title.x = element_text(color="black", size=15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.title.y = element_text(color="black", size=15),
        axis.text.y = element_text(size = 15), 
        legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16))

# ----- Set consistent chromosome order -----
chromosomes <- c("SUPER_1", "SUPER_2", "SUPER_3", "SUPER_X", "SUPER_GRC1", "SUPER_GRC2")

# ----- Plot 1: Total annotated genes per chromosome -----
gtf <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\Annotations\\bcop_core_GRC.gff3",
                comment = "#", col_names = FALSE)
colnames(gtf)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

genes_per_scaffold <- gtf %>%
  filter(feature == "gene") %>%
  mutate(gene_id = str_extract(attribute, "ID=[^;]+") %>% str_remove("ID=")) %>%
  distinct(seqname, gene_id) %>%
  filter(seqname %in% chromosomes) %>%
  count(seqname, name = "Num_genes") %>%
  mutate(Chromosome = fct_relevel(seqname, chromosomes))
genes_per_scaffold

plot_total_genes <- ggplot(genes_per_scaffold, aes(x = Chromosome, y = Num_genes, fill = Chromosome)) +
  geom_col() +
  labs(
    x     = "",
    y     = "Total Annotated Genes",
    title = "A) Total Annotated Genes per Chromosome"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----- Plot 2: Expressed gene counts per chromosome × stage -----
TPM_threshold <- 0.222534
TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\01_STAR_TPM\\combined_TPM_only.tsv", 
                      col_names = FALSE)
colnames(TPM_genes) <- c("Chromosome", "TPM", "Gene", "Species", "Sex", "Tissue", "Stage", "Sample")

counts_per_sample <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% chromosomes) %>%
  distinct(Gene, Stage, Chromosome) %>%
  group_by(Stage, Chromosome) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  mutate(Chromosome = fct_relevel(Chromosome, chromosomes))


plot_expressed <- ggplot(counts_per_sample, aes(x = Chromosome, y = Num_genes, fill = Chromosome)) +
  geom_col() +
  facet_wrap(~ Stage) +
  labs(
    x     = "",
    y     = "Number of Expressed Genes",
    title = "B) Expressed Gene Counts per Chromosome (by Stage)"
  ) +
  box_theme +
  theme(
    legend.position = "nonw",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----- Plot 3: GRC expressed only -----
GRCs  <- c("SUPER_GRC1", "SUPER_GRC2")
GRC_counts_per_sample <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% GRCs) %>%
  distinct(Gene, Stage, Chromosome) %>%
  group_by(Stage, Chromosome) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  mutate(Chromosome = fct_relevel(Chromosome, GRCs))
GRC_counts_per_sample
# Assign consistent colors to each chromosome
chromosome_colors <- scales::hue_pal()(length(chromosomes))
names(chromosome_colors) <- chromosomes

plot_GRC <- ggplot(GRC_counts_per_sample, aes(x = Chromosome, y = Num_genes, fill = Chromosome)) +
  geom_col() +
  facet_wrap(~ Stage) +
  scale_fill_manual(values = chromosome_colors[GRCs]) +
  labs(
    x     = "",
    y     = "Number of Expressed Genes",
    title = "C) Expressed Gene Counts (GRC only)"
  ) +
  box_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_GRC

# ----- Plot 4: Corrected GRC expressed only -----
Corrected_GRC_per_sample <- tibble::tibble(
  Stage = c("0-4h", "0-4h", "4-8h", "4-8h", "adult", "adult", "late-larva-early-pupa", "late-larva-early-pupa"),
  Chromosome = factor(c("SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC1", "SUPER_GRC2",
                        "SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC1", "SUPER_GRC2")),
  Num_genes = c(0, 1, 2, 3, 1, 8, 3, 7)
)
Corrected_GRC_per_sample


plot_GRC_corrected <- ggplot(Corrected_GRC_per_sample, aes(x = Chromosome, y = Num_genes, fill = Chromosome)) +
  geom_col() +
  facet_wrap(~ Stage) +
  scale_fill_manual(values = chromosome_colors[GRCs]) +
  labs(
    x     = "",
    y     = "Number of Expressed Genes",
    title = "D) Expressed GRC gene counts (corrected for mismapping)"
  ) +
  box_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_GRC_corrected

# ----- Plot 4: Expressed vs Total genes per chromosome × stage -----
# Pivot counts_per_sample to wide format to get expression per stage per chromosome
wide_expr <- counts_per_sample %>%
  pivot_wider(names_from = Stage, values_from = Num_genes, values_fill = 0)
wide_expr
# Join total gene counts
combined_counts <- genes_per_scaffold %>%
  rename(Total = Num_genes) %>%
  left_join(wide_expr, by = "Chromosome")

# Reshape to long format with expression stages and total
combined_long <- combined_counts %>%
  select(Chromosome, Total, `0-4h`, `4-8h`, `late-larva-early-pupa`, adult) %>%
  pivot_longer(
    cols = -Chromosome,
    names_to = "Stage",
    values_to = "Count"
  ) %>%
  mutate(
    Stage = factor(Stage, levels = c("Total", "0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, chromosomes)
  )

# Plot
plot_compare <- ggplot(combined_long, aes(x = Chromosome, y = Count, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    x     = "Chromosome",
    y     = "Gene Count",
    title = "A) Total and Expressed Genes per Chromosome"
  ) +
  box_theme +
  scale_fill_manual(
  values = c(
    "Total" = "grey70",
    "0-4h" = "#f3c6f4",      # soft orchid
    "4-8h" = "#e088ec",      # medium orchid
    "late-larva-early-pupa" = "#c043d6",  # deeper violet
    "adult" = "#7a0177"      # rich purple
    ),
    labels = c(
      "Total" = "Total number of genes",
      "0-4h" = "Genes expressed (0–4h)",
      "4-8h" = "Genes expressed (4–8h)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult" = "Genes expressed (adult)"
    ),
    name = NULL  # Optional: remove legend title
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot_compare


# ----- Combine plots -----
expression_plot <- 
  plot_total_genes + plot_expressed +  (plot_GRC +scale_y_continuous(limits =c(0,40))) + (plot_GRC_corrected+scale_y_continuous(limits =c(0,40))) + plot_layout(ncol = 2)


plot_compare / ( (plot_GRC +scale_y_continuous(limits =c(0,40))) +  (plot_GRC_corrected+scale_y_continuous(limits =c(0,40))) )

