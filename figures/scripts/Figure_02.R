library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)

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
GRCs  <- c("SUPER_GRC1", "SUPER_GRC2")

# ----- Plot 1) Gene count v Gene expression per stage -----
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

TPM_threshold <- 0.22

TPM_genes <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\01_STAR_TPM\\combined_TPM_only.tsv", 
                      col_names = FALSE)
colnames(TPM_genes) <- c("Chromosome", "TPM", "Gene", "Species", "Sex", "Tissue", "Stage", "Sample")
TPM_genes

counts_per_sample <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% chromosomes) %>%    # Keep TPM > threshold
  distinct(Gene, Sample, Stage, Chromosome) %>%                   # Keep distinct Gene-Sample pairs (avoid duplicates)
  group_by(Gene, Stage, Chromosome) %>%                           # Group by Gene, Stage, Chromosome
  summarise(Library_count = n(), .groups = "drop") %>%            # Count how many libraries each gene is expressed in
  filter(Library_count >= 2) %>%                                  # Keep genes expressed in at least 2 libraries
  distinct(Gene, Stage, Chromosome) %>%                           # Now get unique Gene-Stage-Chromosome entries
  group_by(Stage, Chromosome) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  mutate(Chromosome = fct_relevel(Chromosome, chromosomes))

counts_per_sample

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
    title = "A"
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
      "0-4h" = "Genes expressed (0–4h embryo)",
      "4-8h" = "Genes expressed (4–8h embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult" = "Genes expressed (adult)"
    ),
    name = NULL  # Optional: remove legend title
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(labels = c(
    "SUPER_1" = "I",
    "SUPER_2" = "II",
    "SUPER_3" = "III",
    "SUPER_X" = "X",
    "SUPER_GRC1" = "GRC1",
    "SUPER_GRC2" = "GRC2"
  )) 

plot_compare


# ----- Plot 2) GRC expression uncorrected -----
# Old way
# GRC_counts_per_sample <- TPM_genes %>%
#   filter(TPM > TPM_threshold, Chromosome %in% GRCs) %>%
#   distinct(Gene, Stage, Chromosome) %>%
#   group_by(Stage, Chromosome) %>%
#   summarise(Num_genes = n(), .groups = "drop") %>%
#   mutate(Chromosome = fct_relevel(Chromosome, GRCs))
# GRC_counts_per_sample

GRC_counts_per_sample <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% GRCs) %>%     # Keep TPM > threshold
  distinct(Gene, Sample, Stage, Chromosome) %>%             # Keep distinct Gene-Sample pairs (avoid duplicates)
  group_by(Gene, Stage, Chromosome) %>%                     # Group by Gene, Stage, Chromosome
  summarise(Library_count = n(), .groups = "drop") %>%      # Count how many libraries each gene is expressed in
  filter(Library_count >= 2) %>%                            # Keep genes expressed in at least 2 libraries
  distinct(Gene, Stage, Chromosome) %>%                     # Now get unique Gene-Stage-Chromosome entries
  group_by(Stage, Chromosome) %>%
  summarise(Num_genes = n(), .groups = "drop") %>%
  mutate(Chromosome = fct_relevel(Chromosome, GRCs))

GRC_counts_per_sample

# Pivot the GRC_counts_per_sample data
GRC_counts_wide <- GRC_counts_per_sample %>%
  pivot_wider(names_from = Stage, values_from = Num_genes, values_fill = 0)

# Reshape to long format
GRC_counts_long <- GRC_counts_wide %>%
  pivot_longer(
    cols = -Chromosome,
    names_to = "Stage",
    values_to = "Count"
  ) %>%
  mutate(
    Stage = factor(Stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, GRCs)
  )

# Plot (uncorrected)
plot_GRC_flat <- ggplot(GRC_counts_long, aes(x = Chromosome, y = Count, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    x = "Chromosome",
    y = "Number of Expressed Genes",
    title = "B"
  ) +
  box_theme +
  scale_fill_manual(
    values = c(
      "0-4h" = "#f3c6f4",
      "4-8h" = "#e088ec",
      "late-larva-early-pupa" = "#c043d6",
      "adult" = "#7a0177"
    ),
    labels = c(
      "0-4h" = "Genes expressed (0–4h Embryo)",
      "4-8h" = "Genes expressed (4–8h Embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult" = "Genes expressed (adult)"
    ),
    name = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +  
  scale_x_discrete(labels = c(
    "SUPER_GRC1" = "GRC1",
    "SUPER_GRC2" = "GRC2"
  )) 

plot_GRC_flat


# ----- Plot 3) GRC expression corrected -----

# Manual tibble of corrected GRC counts
Corrected_GRC_per_sample <- tibble::tibble(
  Stage = c("0-4h", "0-4h", "4-8h", "4-8h", "adult", "adult", "late-larva-early-pupa", "late-larva-early-pupa"),
  Chromosome = factor(c("SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC1", "SUPER_GRC2",
                        "SUPER_GRC1", "SUPER_GRC2", "SUPER_GRC1", "SUPER_GRC2")),
  Num_genes = c(0, 0, 2, 0, 1, 7, 3, 5)
)
Corrected_GRC_per_sample
sum(Corrected_GRC_per_sample$Num_genes)

# Reshape corrected data
Corrected_GRC_long <- Corrected_GRC_per_sample %>%
  mutate(
    Stage = factor(Stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult")),
    Chromosome = fct_relevel(Chromosome, GRCs)
  )

# Plot (corrected)
plot_GRC_corrected <- ggplot(Corrected_GRC_long, aes(x = Chromosome, y = Num_genes, fill = Stage)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    x = "Chromosome",
    y = "Number of Expressed Genes",
    title = "C"
  ) +
  box_theme +
  scale_fill_manual(
    values = c(
      "0-4h" = "#f3c6f4",
      "4-8h" = "#e088ec",
      "late-larva-early-pupa" = "#c043d6",
      "adult" = "#7a0177"
    ),
    labels = c(
      "0-4h" = "Genes expressed (0–4h Embryo)",
      "4-8h" = "Genes expressed (4–8h Embryo)",
      "late-larva-early-pupa" = "Genes expressed (larval/pupa)",
      "adult" = "Genes expressed (adult)"
    ),
    name = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +   
  scale_x_discrete(labels = c(
    "SUPER_GRC1" = "GRC1",
    "SUPER_GRC2" = "GRC2"
  )) 

plot_GRC_corrected

# ----- Combine plots -----

p2 <- plot_compare / ((plot_GRC_flat + scale_y_continuous(limits = c(0, 30))) + 
                  (plot_GRC_corrected + scale_y_continuous(limits = c(0, 30))))

p2
ggsave("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\figures\\combined_bar_plot.svg", 
       plot = p2, width = 17, height = 15, units = "in", dpi = 300)

### Raw numbers ###

#results tables 
genes_per_scaffold
counts_per_sample

# Core genes
core_gene_no <- 6078 + 5020 + 3615 + 2468
core_gene_no 

# GRC
grc_gene_no <- 8101 + 11816
grc_gene_no

# % of genes GRCs comprise
(grc_gene_no / ( grc_gene_no + core_gene_no )) * 100

# % of genes GRCs express
# Number of GRC genes expressed
GRC_expressed <- sum(GRC_counts_per_sample$Num_genes)
# Number of core genes expressed
core_genes_no <- sum((counts_per_sample %>% filter(!Chromosome %in% GRCs))$Num_genes)
GRC_expressed / (GRC_expressed + core_genes_no) * 100


# Total number of genes expressed per stage
core_genes_no
GRC_expressed

# Unique number of genes expressed per stage

GRC_unique_genes <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% GRCs) %>%
  distinct(Gene, Sample, Chromosome) %>%
  group_by(Gene, Chromosome) %>%
  summarise(Library_count = n(), .groups = "drop") %>%
  filter(Library_count >= 2) %>%
  distinct(Gene)

nrow(GRC_unique_genes)

core_unique_genes <- TPM_genes %>%
  filter(TPM > TPM_threshold, Chromosome %in% chromosomes, !Chromosome %in% GRCs) %>%
  distinct(Gene, Sample, Chromosome) %>%
  group_by(Gene, Chromosome) %>%
  summarise(Library_count = n(), .groups = "drop") %>%
  filter(Library_count >= 2) %>%
  distinct(Gene)

nrow(core_unique_genes)



