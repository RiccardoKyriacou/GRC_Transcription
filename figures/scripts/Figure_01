library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggnewscale)
library(patchwork)
library(svglite)

#Theme
box_theme <- theme_bw() + 
  theme(plot.title = element_text(color="black", size=15),
        panel.background = element_rect(fill="white"),
        axis.title.x = element_text(color="black", size=15),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.title.y = element_text(color="black", size=15),
        axis.text.y = element_text(size = 15), 
        legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16)
  )

### Bar plot ###
df <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\03_BLAST_expressed_genes\\GRC_BLAST_table.tsv")

df_counts <- df %>%
  separate(`germ/soma`, into = c("germ_count", "soma_count"), sep = "_", convert = TRUE) %>%
  mutate(
    `%Identity` = replace_na(`%Identity`, 0),
    Coverage = replace_na(Coverage, 0),
    alignment_score = (Coverage * `%Identity`) / 100,
    development_stage = factor(development_stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult"))
  )

# Replace NA in `%Identity` and Coverage with 0
df_counts <- df_counts %>%
  mutate(
    `%Identity` = replace_na(`%Identity`, 0),
    Coverage = replace_na(Coverage, 0),
    alignment_score = (Coverage * `%Identity`) / 100
  )

# Summarize max values per gene and stage
df_summary_counts <- df_counts %>%
  group_by(gene_id, development_stage) %>%
  summarise(
    germ_count = max(germ_count, na.rm = TRUE),
    soma_count = max(soma_count, na.rm = TRUE),
    alignment_score = max(alignment_score, na.rm = TRUE),
    .groups = "drop"
  )

# Reshape for bar plotting
df_long <- df_summary_counts %>%
  pivot_longer(cols = c("germ_count", "soma_count"),
               names_to = "library_type", values_to = "sample_count") %>%
  mutate(sample_type = case_when(
    development_stage %in% c("0-4h", "4-8h") ~ "pre-GRC-elimination",
    library_type == "germ_count" ~ "germline_library",
    library_type == "soma_count" ~ "somatic_library"
  ))


# Set x-position for alignment score color strip
max_sample <- max(df_long$sample_count, na.rm = TRUE)
df_summary_counts <- df_summary_counts %>%
  mutate(tile_x = max_sample + 1)

# Plot
p <- ggplot(df_long, aes(x = sample_count, y = fct_reorder(gene_id, alignment_score), fill = sample_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c(
      "pre-GRC-elimination" = "skyblue",  
      "germline_library" = "#6A5ACD",         
      "somatic_library" = "#F4A261"             
    ),
    name = "Library Type"
  ) +
  labs(
    x = "Sample Count",
    y = "Gene ID",
    title = "Germ vs. Soma Sample Counts per Gene",
    subtitle = "By Developmental Stage"
  ) +
  theme_minimal(base_size = 20) + 
  facet_wrap(~ development_stage, scales = "free_y") 

# Add alignment score color strip
p1 <- p + new_scale_fill() +
  geom_tile(data = df_summary_counts,
            aes(x = tile_x, y = gene_id, fill = alignment_score),
            width = 0.3, height = 0.8, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "inferno", name = "Alignment Score\n(%Coverage × %Identity)")

# Show plot
p1 

ggsave("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\figures\\mismap_plot.svg", 
       plot = p1, width = 25, height = 20, units = "in", dpi = 300)






### Plot with TPM values ##
### Bar plot ###
df <- read_tsv("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\03_BLAST_expressed_genes\\GRC_BLAST_table.tsv")

# Step 1: Parse mean TPMs
df_tpms <- df %>%
  separate(`mean_germ_TPM/mean_soma_TPM`, into = c("mean_germ_TPM", "mean_soma_TPM"), sep = "_", convert = TRUE) %>%
  mutate(
    development_stage = factor(development_stage, levels = c("0-4h", "4-8h", "late-larva-early-pupa", "adult"))
  )

# Step 2: Summarize TPMs (keep highest observed per gene × stage)
df_summary_tpms <- df_tpms %>%
  group_by(gene_id, development_stage) %>%
  summarise(
    mean_germ_TPM = max(mean_germ_TPM, na.rm = TRUE),
    mean_soma_TPM = max(mean_soma_TPM, na.rm = TRUE),
    alignment_score = max((Coverage * `%Identity`) / 100, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Reshape to long format with log TPM
df_tpms_long <- df_summary_tpms %>%
  pivot_longer(cols = c("mean_germ_TPM", "mean_soma_TPM"),
               names_to = "library_type", values_to = "mean_TPM_raw") %>%
  mutate(
    log_TPM = log1p(mean_TPM_raw),
    sample_type = case_when(
      development_stage %in% c("0-4h", "4-8h") ~ "pre-GRC-elimination",
      library_type == "mean_germ_TPM" ~ "germline_library",
      library_type == "mean_soma_TPM" ~ "somatic_library"
    )
  )

# Order genes by max alignment score again
df_tpms_long <- df_tpms_long %>%
  mutate(gene_id = fct_reorder(gene_id, alignment_score))

# Step 4: TPM plot
p_tpm <- ggplot(df_tpms_long, aes(x = log_TPM, y = gene_id, fill = sample_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = c(
      "pre-GRC-elimination" = "#a6c3daff",
      "germline_library" = "#e69f00ff",
      "somatic_library" = "#009e73ff"
    ),
    name = "Library Type"
  ) +
  labs(
    x = "Average Expression log(TPM + 1) per Library Type ",
    y = "Gene ID",
    title = "Germ vs. Soma Mean TPM per Gene (log scale)",
    subtitle = "By Developmental Stage"
  ) +
  theme_minimal(base_size = 20) +
  facet_wrap(~ development_stage, scales = "free_y") 

# Add alignment score 
p_tpm <- p_tpm + new_scale_fill() +
  geom_tile(data = df_summary_counts,
            aes(x = 4.5, y = gene_id, fill = alignment_score),
            width = 0.3, height = 0.8, inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Alignment Score\n(%Coverage × %Identity)")
p_tpm


ggsave("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\figures\\TPM_mismap_plot.svg", 
       plot = p_tpm, width = 25, height = 20, units = "in", dpi = 300)

