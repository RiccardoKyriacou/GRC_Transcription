library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# Data
origin_data <- tribble(
  ~Gene,      ~Category,
  "g11713",   "Transposable Element",
  "g13362",   "Insect",
  "g13363",   "Transposable Element",
  "g13694",   "Insect",
  "g15174",   "Transposable Element",
  "g16029",   "Insect",
  "g17107",   "Transposable Element",
  "g17119",   "Transposable Element",
  "g19161",   "Bacterial",
  "g233",     "Transposable Element",
  "g491",     "Transposable Element",
  "g596",     "Transposable Element",
  "g7610",    "Insect",
  "g7957",    "Insect",
  "g8036",    "Transposable Element"
)

# Count categories and compute label positions
category_counts <- origin_data %>%
  count(Category) %>%
  arrange(n) %>%
  mutate(
    Category = factor(Category, levels = Category),
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    label = paste0(n, "")
  )
category_counts

pie <- ggplot(category_counts, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Category)) +
  geom_rect(color = "white") +
  geom_text(
    aes(
      x = 1.3,
      y = label_pos,
      label = label
    ),
    size = 5,
    hjust = 0.5
  ) +
  coord_polar(theta = "y") +
  xlim(c(-0.5, 1.6)) +
  theme_void() +
  scale_fill_manual(
    labels = c("Bacterial", "Insect", "Transposable Element"),
    values = c(
      "Transposable Element" = "#29335C",
      "Insect" = "#DA70AE",
      "Bacterial" = "#00BFC4"
    )
  ) +
  labs(fill = "Top BLAST Hit") +
  theme(legend.position = "bottom")

pie


# Load the data
df <- read.table(text = "
gene_id tissue stage sex TPM
g11713 germ adult male 0.5
g13362 germ late-larva-early-pupa both-sexes 8.31
g13362 germ adult male 6.61
g13363 germ late-larva-early-pupa male 1.43
g13363 germ adult male 0.9
g13694 germ late-larva-early-pupa male 1.02
g13694 germ adult male 0.66
g15174 germ adult female 0.39
g16029 germ late-larva-early-pupa female 0.57
g17107 germ adult male 1.19
g17119 germ adult female 0.6
g19161 germ late-larva-early-pupa female 1.01
g233 germ late-larva-early-pupa both-sexes 0.45
g491 germ late-larva-early-pupa male 1.06
g596 germ late-larva-early-pupa female 1.26
g7610 germ adult both-sexes 0.26
g7957 germ 4-8h female 0.44
g8036 germ 4-8h female 0.34
", header = TRUE)

# Pivot data for heatmap
df_heat <- df %>%
  group_by(gene_id, stage) %>%
  summarise(TPM = sum(TPM), .groups = "drop") %>%
  pivot_wider(names_from = stage, values_from = TPM, values_fill = 0)

# Convert to long format for ggplot
df_long <- df_heat %>%
  pivot_longer(-gene_id, names_to = "stage", values_to = "TPM")


# Join category info to df_long
df_long_ordered <- df_long %>%
  left_join(origin_data, by = c("gene_id" = "Gene")) %>%
  mutate(
    stage = factor(stage, levels = c("4-8h", "late-larva-early-pupa", "adult")),
    Category = factor(Category, levels = c("None", "Bacterial", "Insect", "Transposable Element"))
  ) %>%
  arrange(Category, gene_id) %>%
  mutate(gene_id = factor(gene_id, levels = rev(unique(gene_id))))


heat <- ggplot(df_long_ordered, aes(x = stage, y = gene_id, fill = log2(TPM + 1))) +
  geom_tile(color = "grey30") +
  scale_fill_gradientn(
    colors = c("white", "#F8C464", "#FC8D59", "#E34A33", "#D93838", "#840715"),
    name = "log2(TPM + 1)"
  )+
  theme_minimal(base_size = 12) +
  labs(title = "TPM Expression Heatmap", x = "Life Stage", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heat
heat <- ggplot(df_long_ordered, aes(x = stage, y = gene_id, fill = log2(TPM+1))) +
  geom_tile(color = "grey30") +
  scale_fill_gradientn(
    colors = c("white", "#7F66E0", "#4B30C2"),
    name = "log2(TPM + 1)"
  )+
  theme_minimal(base_size = 12) +
  labs(title = "TPM Expression Heatmap", x = "Life Stage", y = "Gene") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heat
pie_heat <- pie + heat 
pie_heat

ggsave("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\GRC_expression_analysis\\B_coprophila\\figures\\new_heatmap_pie_plot.svg", 
       plot = pie_heat, width = 15, height = 15, units = "in", dpi = 300)

#TODO change heatmap order

