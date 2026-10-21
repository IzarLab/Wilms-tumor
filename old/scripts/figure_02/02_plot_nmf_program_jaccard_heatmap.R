# ============================================================
# Figure 2 NMF metaprogram Jaccard heatmap
# ============================================================
#
# GitHub script name : 02_plot_nmf_program_jaccard_heatmap.R
# Original file      : Fig.2.NMF_program_heatmap_cleaned(1).R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Uses NMF program Jaccard matrix and manually curated program annotations to generate cleaned metaprogram heatmap.
#
# Usage:
#   Run from the root of the GitHub repository, or edit the input/output
#   paths in the SETTINGS / Input-output section below before running.
#
# Notes for reproducibility:
#   - This script preserves the analytical logic of the original working
#     manuscript script.
#   - Hard-coded local paths from the original analysis may need to be
#     changed to repository-relative paths such as data/processed/ and
#     results/figure_XX/.
#   - Large objects such as Seurat .RData files should usually be stored
#     outside GitHub or tracked with Git LFS / external data release.
#
# Suggested output folder:
#   results/figure_02/
#
# ============================================================

# This script plots a cleaned Jaccard similarity matrix of clustered factors where each cluster is called a program.

# Input: It needs Jaccard_similarity_matrix_programs.tsv and manually annotated programs_annotation_processed.tsv from NMF_metaprogram_heatmap.R.
# Output: It outputs a cleaned heatmap (NMF_programs_cleaned.pdf) of clustered factors where each cluster is called a program.
#         The file programs_annotation.tsv MUST BE ANNOTATED (now named programs_annotation_processed.tsv) USIGING THE OUTPUT
#         HEATMAP (programs.pdf) AND NMF_GSEA_1e-3.tsv.

# NMF-specific terminology:

# Feature: genes
# Sample: cells (in single cell data)
# Factor: basis component (similar to principle component) or a metagene
# Rank: max number of factors (similar to rank in prcomp() for PCA)
# V: data matrix: genes × cells (unlike prcomp that input is cells x genes)
# W: basis matrix (similar to rotation/loadings matrix in PCA): genes × factors
# H: Coefficient Matrix (similar to x matrix in PCA but always non-negative): factors × samples (cells). Unlike PCA, each value is a metagene expression not a coordinate!

library(pheatmap)
library(gridExtra)
library(viridis)
set.seed(42)

# Reading in data ####

# reading in raw Jaccard similarity matrix
mat_ = as.matrix(read.delim(file = '~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2/Jaccard_similarity_matrix_programs.tsv', header = T, sep = '\t', as.is = T, check.names = F))

# reading in modified annotation file
# Att.: FIRST ANNOTATE open programs_annotation.tsv using NMF_programs.pdf and NMF_GSEA_1e-3.tsv files

annot_ = read.delim(file = '~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2/programs_annotation_processed.tsv', header = T, sep = '\t', as.is = T, check.names = T)     # manually annotated programs_annotation.tsv

# programs identified in programs_annotation_processed.tsv for ordering levels
progs_ = c('blastema 1','S blastema 1','G2M blastema 1',
           'blastema 2',
           'PAX3+ myogenic precursors','striated myocyte',
           'smooth-myocyte precursors','nascent smooth myocyte','smooth myocyte',
           'podocyte',
           'tubular','connecting tubule')

# reformatting annotation file for ggpolt2
ls_ = list()
for(p_ in progs_)
{
  ls_[[p_]] = annot_[annot_$program %in% p_,]
  ls_[[p_]] = ls_[[p_]][order(ls_[[p_]]$condition),]
}
annot_ = do.call(ls_, what = rbind)
rownames(annot_) = annot_$factor_sample
annot_$program = factor(x = annot_$program, levels = progs_)

# plotting final Jaccard heatmap

pdf(file = '~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2/NMF_programs_cleaned.pdf', width = 10, height = 10)

cols_1 = sample(viridis(n = length(progs_), option = 'turbo'))                    # color for each program
names(cols_1) = levels(annot_$program)

cols_2 = sample(viridis(n = length(unique(annot_$sample)), option = 'turbo'))     # color for each patient sample
names(cols_2) = unique(annot_$sample)

cols_3 = c(untreated = 'blue', neoadjuvant = 'red')                               # color for treatment condition
annot_cols = list(program = cols_1, sample = cols_2, condition = cols_3)

m_ = mat_[rownames(annot_), rownames(annot_)]
m_[m_ <= 0.05] = 0
p_ = pheatmap(mat = m_, annotation_row = annot_[,-1, drop = F], annotation_col = annot_[,-1, drop = F],
              cluster_cols = F, cluster_rows = F, clustering_method = 'ward.D2', treeheight_row = 0, treeheight_col = 0,
              border_color = NA, scale = 'none', color = rev(viridis(n = 100, option = 'A')),
              show_rownames = F, show_colnames = F, fontsize = 2,
              annotation_colors = annot_cols,
              legend_breaks = c(0,.5,1),
              silent = F)
p_
graphics.off()


##### NMF
W944 <- load('~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2/NMF_W944.RData')
W944 <- r_estimates
colnames(W944$measures)
table(W944$measures$rank)

names(W944)
class(W944)

library(NMF)

##### output directory
outdir <- "W944_NMF_WH_matrices"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

##### available ranks
ranks <- names(W944$fit)

##### loop through all ranks
for (k in ranks) {
  
  cat("Processing rank:", k, "\n")
  
  fit_k <- W944$fit[[k]]
  
  ##### try direct extraction
  W <- tryCatch(
    basis(fit_k),
    error = function(e) NULL
  )
  
  H <- tryCatch(
    coef(fit_k),
    error = function(e) NULL
  )
  
  ##### fallback: embedded fit slot
  if (is.null(W) || is.null(H) || is.null(dim(W)) || is.null(dim(H))) {
    
    cat("Using embedded @fit slot for rank", k, "\n")
    
    model_k <- fit_k@fit
    
    W <- basis(model_k)
    H <- coef(model_k)
  }
  
  ##### dimensions
  cat("W dimensions:", dim(W), "\n")
  cat("H dimensions:", dim(H), "\n")
  
  ##### save csv
  write.csv(
    W,
    file = file.path(
      outdir,
      paste0("rank_", k, "_W_genes_by_programs.csv")
    )
  )
  
  write.csv(
    H,
    file = file.path(
      outdir,
      paste0("rank_", k, "_H_programs_by_cells.csv")
    )
  )
  
  ##### save RDS
  saveRDS(
    W,
    file = file.path(
      outdir,
      paste0("rank_", k, "_W.rds")
    )
  )
  
  saveRDS(
    H,
    file = file.path(
      outdir,
      paste0("rank_", k, "_H.rds")
    )
  )
}

cat("Done.\n")

##### NMF annotations - downstream
#Integrate factors

library("dplyr")
library("plyr")
library("readr")
library("purrr")
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(cluster)

data_all <- list.files(path = "./nmf_matrices/", pattern = "*.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "genes")
data_all[is.na(data_all)] <- 0

write.csv(data_all,file="~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2/wilms_tumor_gene_top100_correlation.csv")

############################################################
# NMF W-matrix top unique marker dotplot
#
# Input:
#   nmf-dotplot.csv
#
# Assumption:
#   - First column = gene names
#   - Remaining columns = annotated metaprograms
#   - Values = NMF W weights
#
# Output:
#   - Continuous y-axis dotplot
#   - Exactly top_n markers per metaprogram where possible
#   - Unique markers assigned to only one metaprogram
#   - PDF, SVG, and marker CSV
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

############################################################
# User parameters
############################################################

infile <- "nmf-dotplot.csv"
outdir <- "NMF_dotplot_unique_continuous"
top_n <- 5   # change to 5 if needed

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
# Read W matrix
############################################################

w_raw <- read.csv(
  infile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# First column is assumed to be gene name
colnames(w_raw)[1] <- "gene"

w_raw <- w_raw %>%
  dplyr::mutate(gene = as.character(gene)) %>%
  dplyr::filter(!is.na(gene), gene != "")

annotation_order <- colnames(w_raw)[-1]

############################################################
# Convert W matrix to long format
############################################################

w_long <- w_raw %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "annotation",
    values_to = "weight"
  ) %>%
  dplyr::mutate(
    annotation = as.character(annotation),
    gene = as.character(gene),
    weight = as.numeric(weight)
  ) %>%
  dplyr::filter(
    !is.na(annotation),
    !is.na(gene),
    !is.na(weight)
  )

############################################################
# Rank genes within each metaprogram
############################################################

w_ranked <- w_long %>%
  dplyr::mutate(
    annotation_factor = factor(annotation, levels = annotation_order)
  ) %>%
  dplyr::arrange(
    annotation_factor,
    dplyr::desc(weight)
  )

############################################################
# Select unique top markers per metaprogram
#
# This assigns each gene to only one metaprogram.
# Earlier metaprograms in annotation_order get priority.
# Increase candidate_pool_factor if some groups have < top_n markers.
############################################################

candidate_pool_factor <- 5

selected_list <- list()
used_genes <- character(0)

for (ann in annotation_order) {
  
  tmp <- w_ranked %>%
    dplyr::filter(annotation == ann) %>%
    dplyr::arrange(dplyr::desc(weight)) %>%
    dplyr::filter(!gene %in% used_genes) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      rank_in_annotation = dplyr::row_number(),
      block_annotation = ann
    )
  
  selected_list[[ann]] <- tmp
  used_genes <- unique(c(used_genes, tmp$gene))
}

top_markers <- dplyr::bind_rows(selected_list) %>%
  dplyr::mutate(
    annotation = factor(annotation, levels = annotation_order),
    block_annotation = factor(block_annotation, levels = annotation_order)
  ) %>%
  dplyr::arrange(annotation, rank_in_annotation)

############################################################
# Check marker counts per metaprogram
############################################################

marker_counts <- top_markers %>%
  dplyr::count(annotation, name = "n_markers")

print(marker_counts)

if (any(marker_counts$n_markers < top_n)) {
  warning(
    "Some metaprograms have fewer than top_n unique markers. ",
    "This happens when high-ranking genes overlap heavily across metaprograms."
  )
}

############################################################
# Continuous block-wise y-axis ordering
############################################################

gene_order <- top_markers %>%
  dplyr::arrange(annotation, rank_in_annotation) %>%
  dplyr::pull(gene)

top_markers <- top_markers %>%
  dplyr::mutate(
    gene_block = factor(gene, levels = rev(gene_order))
  )

############################################################
# Save marker table
############################################################

write.csv(
  top_markers,
  file = file.path(
    outdir,
    paste0("NMF_top", top_n, "_unique_markers_per_metaprogram.csv")
  ),
  row.names = FALSE
)

write.csv(
  marker_counts,
  file = file.path(
    outdir,
    paste0("NMF_top", top_n, "_marker_counts_per_metaprogram.csv")
  ),
  row.names = FALSE
)

############################################################
# Plot: continuous y-axis, grey-to-black gradient
############################################################

p <- ggplot(
  top_markers,
  aes(
    x = annotation,
    y = gene_block,
    size = weight,
    color = weight
  )
) +
  geom_point() +
  scale_color_gradient(
    low = "grey85",
    high = "black",
    name = "NMF W\nweight"
  ) +
  scale_size_continuous(
    range = c(1.5, 5),
    name = "NMF W\nweight"
  ) +
  labs(
    x = "Annotated metaprograms",
    y = "Unique marker genes ordered by metaprogram block",
    title = paste0("Top ", top_n, " unique markers per annotated metaprogram")
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 10),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    ),
    legend.position = "right"
  )

############################################################
# Save plots
############################################################

plot_height <- max(5, 0.18 * nrow(top_markers))

ggsave(
  filename = file.path(
    outdir,
    paste0("NMF_top", top_n, "_unique_marker_dotplot_continuous.pdf")
  ),
  plot = p,
  width = 4,
  height = plot_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    outdir,
    paste0("NMF_top", top_n, "_unique_marker_dotplot_continuous.svg")
  ),
  plot = p,
  width = 4,
  height = plot_height,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

############################################################
# Optional biological sanity table:
# Top genes collapsed per annotation
############################################################

bio_summary <- top_markers %>%
  dplyr::group_by(annotation) %>%
  dplyr::summarise(
    top_genes = paste(gene, collapse = ", "),
    .groups = "drop"
  )

write.csv(
  bio_summary,
  file = file.path(
    outdir,
    paste0("NMF_top", top_n, "_genes_by_metaprogram_summary.csv")
  ),
  row.names = FALSE
)

cat("Done.\n")
cat("Outputs saved in:", outdir, "\n")