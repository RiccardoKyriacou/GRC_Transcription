# Plot chromosome alignments as dot plots
################################################################################
require(dplyr)
require(ggplot2)
require(S4Vectors)
require(stringr)
library(patchwork)
################################################################################

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#iocManager::install("S4Vectors")

# Functions
# read paf alignment file
load_paf_file <- function(file, t_chroms = NULL, t_ranges = NULL, t_sp = NULL,
                          q_chroms = NULL, q_ranges = NULL, q_sp = NULL){
  
  # read file
  ali <- read.table(file = file, header = FALSE)[,c(1:9)]
  colnames(ali) <- c("target", "t_length", "t_start", "t_end", "t_string",
                     "query", "q_length", "q_start", "q_end")
  
  # filter target chromosomes
  if(!is.null(t_chroms)){
    ali <- ali %>%
      filter(target %in% t_chroms)
  }
  
  # filter query chromosomes
  if(!is.null(q_chroms)){
    ali <- ali %>%
      filter(query %in% q_chroms)
  }
  
  # filter target ranges
  if(!is.null(t_ranges)){
    ali <- ali %>%
      filter((t_start > t_ranges[1]) & (t_end < t_ranges[2]))
  }
  
  # filter q_ranges
  if(!is.null(q_ranges)){
    ali <- ali %>%
      filter((q_start > q_ranges[1]) & (q_end < q_ranges[2]))
  }
  
  ali <- ali %>% arrange(target, t_start)
  ali$t_species <- rep(t_sp, nrow(ali))
  ali$q_species <- rep(q_sp, nrow(ali))
  return(ali)
}

# read chromosome files
load_chrom_sizes <- function(file, species){
  df <- read.csv(file, header = FALSE, sep = "\t")
  colnames(df) <- c("chrom", "size")
  df["species"] <- rep(species, length(df$chrom))
  return(df)
}

# convert genomic coordinates to linear format
prepare_linear_coordinates <- function(ali, t_sp, q_sp, chrom_list){
  # read chromosome sizes
  chrom_sizes <- list()
  for(sp in c(t_sp, q_sp)){
    chrom_sizes[[sp]] <- load_chrom_sizes(file = chrom_list[[sp]], species = sp)
  }
  
  # merge chromosomes sizes with the alignment table
  ali <- left_join(ali, chrom_sizes[[t_sp]], 
                   by = c("target" = "chrom", "t_species" = "species"))
  
  ali <- left_join(ali, chrom_sizes[[q_sp]],
                   by = c("query" = "chrom", "q_species" = "species"))
  
  names(ali)[names(ali) == "size.x"] <- "t_chrom_size"
  names(ali)[names(ali) == "size.y"] <- "q_chrom_size"
  
  # transform coordinates of the target to linear format
  ali <- ali %>% arrange(target, t_start)
  
  chr_offset <- 0
  ali_lin <- NULL
  
  for(i in unique(ali$target)){
    ali_t <- ali[ali$target == i,]
    size_t <- unique(ali_t$t_chrom_size)
    
    ali_t$t_linear_start <- chr_offset + ali_t$t_start
    ali_t$t_linear_end <- chr_offset + ali_t$t_end
    ali_lin <- rbind(ali_lin, ali_t)
    
    chr_offset <- chr_offset + size_t
  }
  
  # transform coordinates of the query to linear format
  ali_lin <- ali_lin %>% arrange(query, q_start)
  
  chr_offset <- 0
  ali_lin_both <- NULL
  
  for(i in unique(ali_lin$query)){
    ali_q <- ali_lin[ali_lin$query == i,]
    size_q <- unique(ali_q$q_chrom_size)
    
    ali_q$q_linear_start <- chr_offset + ali_q$q_start
    ali_q$q_linear_end <- chr_offset + ali_q$q_end
    ali_lin_both <- rbind(ali_lin_both, ali_q)
    
    chr_offset <- chr_offset + size_q
  }
  
  ali_lin_both <- ali_lin_both %>% arrange(target, t_linear_start) %>%
    select(target, t_linear_start, t_linear_end, t_species, t_chrom_size,
           query, q_linear_start, q_linear_end, q_species, q_chrom_size)
  
  return(ali_lin_both)
}

# calculate break positions for plotting
zipWithNext <- function(x, step = 1) {
  Pairs(
    x,
    c(tail(x, step * -1), rep(NA, step))
  )
}

calcLabelPosition <- function(breakPos) {
  # Return the position of the midpoint of each chr
  # in the context of the merged object
  breakPos |> zipWithNext() |> as.data.frame() |>
    rowMeans() |> head(-1)
}
################################################################################
# Set up main info here

setwd("C:\\Users\\s2673271\\OneDrive - University of Edinburgh\\PhD\\Y1\\Sciaridae\\long_reads\\FastGA")
home <- getwd()
home

target_species <- "Bcop"
# target chromosomes to plot
target_chrom2plot <- c("651kb_full_HGT_region")


query_species <- "Rickettsia"
# query chromosomes to plot
query_chrom2plot <- c("ptg000159c")

chrom_files <- c("651kb_length.tsv",
                 "Rickettsiaceae_contig_sizes.tsv")

chrom_list <- file.path(home, chrom_files)
names(chrom_list) <- c("Bcop", "Rickettsia")

alignment_file <- file.path(home, "grc2_HGTregion_vs_rickettsia_1to1.1aln.paf")

################################################################################
### --- load paf alignment --- ###
ali <- load_paf_file(
  file = alignment_file,
  t_chroms = target_chrom2plot, t_sp = target_species, 
  q_chroms = query_chrom2plot, q_sp = query_species)

### --- convert coordinates to linear format --- ###
ali_lin <- prepare_linear_coordinates(
  ali = ali, t_sp = target_species, q_sp = query_species, chrom_list = chrom_list)

### --- plotting --- ###
p <- ggplot(ali_lin)  +
  aes(x = t_linear_start, y = q_linear_start, xend = t_linear_end, yend = q_linear_end)

p <- p + geom_segment(aes(colour = target), 
                      lineend = "round", linewidth = 1.2, colour="#00BFC4") + 
  labs(x = target_species, y = query_species) +
  theme_bw()
p

### Plot wide ####
# calculate breaks and labels positions
breaks <- list()
ticks  <- list()
labels <- list()

chr_info_t <- ali_lin %>%
  arrange(target) %>%
  distinct(target, t_chrom_size) %>%
  mutate(cum_end = cumsum(t_chrom_size))

chr_info_q <- ali_lin %>%
  arrange(query) %>%
  distinct(query, q_chrom_size) %>%
  mutate(cum_end = cumsum(q_chrom_size))

breaks$target <- c(0, chr_info_t$cum_end)
ticks$target <- calcLabelPosition(breaks$target)
labels$target <- chr_info_t$target

# add breaks and labels
# total size of target chromosome(s)
target_total <- max(chr_info_t$cum_end)

p_wide <- p +
  scale_x_continuous(expand = c(0, 0), minor_breaks = NULL,
                     limits = c(0, target_total),
                     breaks = breaks$target, labels = NULL, position = 'top',
                     sec.axis = dup_axis(breaks=ticks$target, labels=labels$target)) +
  scale_y_continuous(expand = c(0, 0), minor_breaks = NULL,
                     breaks = breaks$query, labels = NULL, position = 'right',
                     sec.axis = dup_axis(breaks=ticks$query, labels=labels$query)) +
  coord_fixed(ratio = 1/10, clip = "off") +   # squish vertically
  theme(legend.position="none")

p_wide

ggsave(plot = p_wide, filename = file.path(home, "651kb_HGT_region_dotplot.svg"), 
       device = "svg")#, units = "cm", width = 100, height = 40)

# compute total ranges from your chr_info objects
#target_total <- max(chr_info_t$cum_end)   # 651000 in your example
#query_total  <- max(chr_info_q$cum_end)   # 1897780 in your example

# ratio to pass to coord_fixed
#ratio <- target_total / query_total

#p_wide <- p +
#  scale_x_continuous(expand = c(0, 0), limits = c(0, target_total),
#                     breaks = breaks$target, labels = NULL, position = 'top',
#                     sec.axis = dup_axis(breaks = ticks$target, labels = labels$target)) +
#  scale_y_continuous(expand = c(0, 0),
#                     breaks = breaks$query, labels = NULL, position = 'right',
#                     sec.axis = dup_axis(breaks = ticks$query, labels = labels$query)) +
#  coord_fixed(ratio = ratio, clip = "off") +
#  theme(legend.position = "none")


################################################################################
#Pie chart showing % bacteria gene integrated into GRC2
################################################################################

library(dplyr)
library(ggplot2)

# ---- Input ----
paf_file <- "grc2_HGTregion_vs_rickettsia_1to1.1aln.paf"
query_genome_size <- 1753321 + 98675 # primary contigs only = ptg000159c + ptg000322c

# ---- Load PAF ----
paf <- read.table(paf_file, header = FALSE, stringsAsFactors = FALSE)
colnames(paf)[1:9] <- c("target","t_len","t_start","t_end","strand",
                        "query","q_len","q_start","q_end")

# ---- Calculate total aligned query bases ----
aligned_bases <- sum(paf$q_end - paf$q_start + 1)
aligned_bases
perc <- 100 * aligned_bases / query_genome_size
cat("Aligned bases:", aligned_bases, "\n")
cat("Genome size:", query_genome_size, "\n")
cat("Percentage integrated:", round(perc, 2), "%\n")

# ---- Pie chart ----
# ---- Prep for pie chart ----
df_pie <- data.frame(
  Category = c("Integrated", "Not integrated"),
  bases = c(aligned_bases, query_genome_size - aligned_bases)
)

# cumulative sums for ymin/ymax
df_pie <- df_pie %>%
  arrange(desc(Category)) %>%
  mutate(
    fraction = bases / sum(bases),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    label_pos = (ymax + ymin) / 2,
    label = paste0(round(fraction*100, 1), "%")
  )

# ---- Plot ----
pie <- ggplot(df_pie, aes(ymax=ymax, ymin=ymin, xmax=1, xmin=0, fill=Category)) +
  geom_rect(color="white") +
  geom_text(
    aes(
      x=1.3,
      y=label_pos,
      label=label
    ),
    size=10,
    hjust=0.5
  ) +
  coord_polar(theta="y") +
  xlim(c(-0.5, 1.6)) +
  theme_void() +
  scale_fill_manual(
    values = c(
      "Integrated" = "#00BFC4",
      "Not integrated" = ""
    )
  ) +
  labs(
    fill = "Rickettsia genome integrated into GRC2"
  ) +
  theme(legend.position="right")

pie

ggsave(plot = pie, filename = file.path(home, "HGT_piechart.svg"), 
       device = "svg")#, units = "cm", width = 100, height = 40)
