#!/usr/bin/env Rscript

# =============================================================================
# Script: 02_plot_global_umaps.R
# Project: Wilms tumor plasticity manuscript
# Figure: Figure 1 - global TME UMAP panels
#
# Purpose:
#   Generate UMAP plots from the full Wilms tumor TME Seurat object, colored by
#   major compartment, treatment condition, and cell-cycle phase.
#
# Input:
#   data/processed/TME.RData
#       Seurat object named `tme_` with UMAP coordinates and metadata columns:
#       compartment_name, condition, Phase.
#
# Outputs:
#   figures/figure_01/global_umap_compartment_no_label.pdf
#   figures/figure_01/global_umap_compartment_label.pdf
#   figures/figure_01/global_umap_condition_no_label.pdf
#   figures/figure_01/global_umap_condition_label.pdf
#   figures/figure_01/global_umap_phase_no_label.pdf
#   figures/figure_01/global_umap_phase_label.pdf
#
# Run from repository root:
#   Rscript scripts/figure_01/02_plot_global_umaps.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# User-configurable paths
# -----------------------------------------------------------------------------
tme_path <- file.path("data", "processed", "TME.RData")
out_dir <- file.path("figures", "figure_01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tme_path)) {
  stop("Missing required input file: ", tme_path, call. = FALSE)
}

load(tme_path)
if (!exists("tme_")) {
  stop("The TME.RData file must contain a Seurat object named `tme_`.", call. = FALSE)
}

obj <- tme_

required_cols <- c("compartment_name", "condition", "Phase")
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing columns in obj@meta.data: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# Factor ordering and color palettes
# -----------------------------------------------------------------------------
obj$compartment_name <- factor(obj$compartment_name, levels = c("cancer", "immune", "vasculature"))
obj$condition <- factor(obj$condition, levels = c("untreated", "neoadjuvant"))
obj$Phase <- factor(obj$Phase, levels = c("G1", "S", "G2M"))

compartment_cols <- c(
  "cancer" = "#D73027",
  "immune" = "#4575B4",
  "vasculature" = "#1B9E77"
)

condition_cols <- c(
  "untreated" = "#7F7F7F",
  "neoadjuvant" = "#E69F00"
)

phase_cols <- c(
  "G1" = "#66C2A5",
  "S" = "#FC8D62",
  "G2M" = "#8DA0CB"
)

# -----------------------------------------------------------------------------
# Plotting helper
# -----------------------------------------------------------------------------
plot_umap <- function(obj, group_by, cols, title, label = FALSE, pt_size = 0.1) {
  DimPlot(
    obj,
    reduction = "umap",
    group.by = group_by,
    label = label,
    repel = label,
    label.size = 5,
    pt.size = pt_size,
    raster = FALSE
  ) +
    scale_color_manual(values = cols, drop = FALSE) +
    ggtitle(title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

# -----------------------------------------------------------------------------
# Generate and save plots
# -----------------------------------------------------------------------------
plots <- list(
  global_umap_compartment_no_label = plot_umap(obj, "compartment_name", compartment_cols, "UMAP - Compartment", label = FALSE),
  global_umap_compartment_label = plot_umap(obj, "compartment_name", compartment_cols, "UMAP - Compartment", label = TRUE),
  global_umap_condition_no_label = plot_umap(obj, "condition", condition_cols, "UMAP - Condition", label = FALSE),
  global_umap_condition_label = plot_umap(obj, "condition", condition_cols, "UMAP - Condition", label = TRUE),
  global_umap_phase_no_label = plot_umap(obj, "Phase", phase_cols, "UMAP - Cell-cycle phase", label = FALSE),
  global_umap_phase_label = plot_umap(obj, "Phase", phase_cols, "UMAP - Cell-cycle phase", label = TRUE)
)

for (plot_name in names(plots)) {
  out_file <- file.path(out_dir, paste0(plot_name, ".pdf"))
  ggsave(out_file, plots[[plot_name]], width = 8, height = 6)
  message("Saved: ", out_file)
}
