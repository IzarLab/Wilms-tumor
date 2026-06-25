# ============================================================
# Figure 5 immune/vascular subtype marker dot plots
# ============================================================
#
# GitHub script name : 02_plot_tme_subtype_marker_dotplots.R
# Original file      : dotplots.wims.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Computes top markers and generates dot plots for lymphoid, myeloid, and vasculature subclusters.
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
#   results/figure_05/
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# 1. Load objects
# =========================================================
load("~/Documents/Izar_Group/wilms/integrated_compartment_lymphoid_processed.RData")
lymphoid <- s_objs

load("~/Documents/Izar_Group/wilms/integrated_compartment_vasculature_processed.RData")
vasculature <- s_objs

load("~/Documents/Izar_Group/wilms/integrated_compartment_myeloid_processed.RData")
myeloid <- s_objs

# =========================================================
# 2. Output directory
# =========================================================
outdir <- "~/Documents/Izar_Group/wilms/dotplots_top5_markers"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 3. Helper function
# =========================================================
compute_and_save_markers <- function(
    obj,
    group_col = "cell_type",
    assay_use = NULL,
    slot_use = "data",
    min_pct = 0.1,
    logfc_thresh = 0.25,
    test_use = "wilcox",
    only_pos = TRUE,
    out_rds = "markers.rds",
    out_csv = "markers.csv"
) {
  stopifnot(group_col %in% colnames(obj@meta.data))
  
  if (is.null(assay_use)) assay_use <- DefaultAssay(obj)
  if ("RNA" %in% Assays(obj)) assay_use <- "RNA"
  DefaultAssay(obj) <- assay_use
  
  obj[[group_col]] <- as.character(obj[[group_col, drop = TRUE]])
  Idents(obj) <- obj[[group_col, drop = TRUE]]
  
  markers <- FindAllMarkers(
    object = obj,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_thresh,
    test.use = test_use,
    assay = assay_use,
    slot = slot_use
  )
  
  saveRDS(markers, out_rds)
  write.csv(markers, out_csv, row.names = FALSE)
  
  invisible(markers)
}

# ggsave(
#     filename = pdf_file,
#     plot = p,
#     width = width,
#     height = height,
#     units = "in"
#   )
#   
#   # also save marker table
#   out_csv <- sub("\\.pdf$", "_top_markers.csv", pdf_file)
#   write.csv(top_markers, out_csv, row.names = FALSE)
#   
#   invisible(list(plot = p, markers = markers, top_markers = top_markers))
# }

`%||%` <- function(a, b) if (!is.null(a)) a else b

#--------------------------------------------------
# Helper: plot from precomputed markers
#--------------------------------------------------
plot_dotplot_from_markers <- function(
    obj,
    markers,
    group_col = "cell_type",
    group_order = NULL,
    top_n = 5,
    title = "Top markers per subgroup",
    pdf_file = "dotplot.pdf",
    assay_use = NULL,
    dot_scale = 6,
    low_col = "grey90",
    high_col = "black",
    width = 8,
    height = 14
) {
  stopifnot(group_col %in% colnames(obj@meta.data))
  
  if (is.null(assay_use)) assay_use <- DefaultAssay(obj)
  if ("RNA" %in% Assays(obj)) assay_use <- "RNA"
  DefaultAssay(obj) <- assay_use
  
  obj[[group_col]] <- as.character(obj[[group_col, drop = TRUE]])
  
  current_groups <- sort(unique(obj[[group_col, drop = TRUE]]))
  if (is.null(group_order)) {
    group_order <- current_groups
  } else {
    group_order <- intersect(group_order, current_groups)
  }
  
  obj[[group_col]] <- factor(obj[[group_col, drop = TRUE]], levels = group_order)
  Idents(obj) <- obj[[group_col, drop = TRUE]]
  
  fc_col <- c("avg_log2FC", "avg_logFC")[c("avg_log2FC", "avg_logFC") %in% colnames(markers)][1]
  if (is.na(fc_col)) stop("Could not find avg_log2FC / avg_logFC in marker table.")
  
  # keep only requested groups
  markers2 <- markers %>%
    filter(cluster %in% group_order)
  
  # top N markers per subgroup
  top_markers <- markers2 %>%
    group_by(cluster) %>%
    arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # preserve block order
  top_markers$cluster <- factor(top_markers$cluster, levels = group_order)
  top_markers <- top_markers %>%
    arrange(cluster, desc(.data[[fc_col]]))
  
  # unique genes while keeping first occurrence
  gene_tbl <- top_markers %>%
    select(cluster, gene, all_of(fc_col)) %>%
    distinct() %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup()
  
  # make explicit block structure with spacer rows
  gene_blocks <- lapply(seq_along(group_order), function(i) {
    g <- group_order[i]
    genes_g <- gene_tbl %>%
      filter(cluster == g) %>%
      pull(gene)
    
    if (length(genes_g) == 0) return(NULL)
    
    out <- data.frame(
      feature = genes_g,
      block = g,
      is_spacer = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (i < length(group_order)) {
      out <- rbind(
        out,
        data.frame(
          feature = paste0("___SPACER___", g),
          block = g,
          is_spacer = TRUE,
          stringsAsFactors = FALSE
        )
      )
    }
    out
  })
  
  gene_layout <- bind_rows(gene_blocks)
  feature_order <- rev(gene_layout$feature)
  
  real_features <- gene_layout$feature[!gene_layout$is_spacer]
  
  # build dotplot only on real genes
  dp <- DotPlot(
    object = obj,
    features = rev(real_features),
    group.by = group_col,
    assay = assay_use,
    cols = c(low_col, high_col),
    dot.scale = dot_scale
  )
  
  plot_df <- dp$data
  
  # rename for clarity
  plot_df$cell_group <- factor(plot_df$id, levels = group_order)
  plot_df$feature <- as.character(plot_df$features.plot)
  
  # insert spacer rows manually
  full_y_levels <- feature_order
  plot_df$feature <- factor(plot_df$feature, levels = full_y_levels)
  
  # subgroup labels placed near spacer lines
  block_centers <- gene_layout %>%
    mutate(y_id = rev(seq_len(n()))) %>%
    filter(!is_spacer) %>%
    group_by(block) %>%
    summarise(y = mean(y_id), .groups = "drop")
  
  p <- ggplot(plot_df, aes(x = cell_group, y = feature)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_size(range = c(0, dot_scale)) +
    scale_color_gradient(low = low_col, high = high_col) +
    theme_bw(base_size = 11) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.margin = margin(10, 80, 10, 10)
    ) +
    ggtitle(title) +
    coord_cartesian(clip = "off")
  
  # horizontal separators after each block
  spacer_pos <- gene_layout %>%
    mutate(y_id = rev(seq_len(n()))) %>%
    filter(is_spacer) %>%
    pull(y_id)
  
  if (length(spacer_pos) > 0) {
    p <- p + geom_hline(
      yintercept = spacer_pos + 0.5,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.4
    )
  }
  
  # right-side block labels
  label_df <- block_centers
  label_df$x <- length(group_order) + 0.6
  
  p <- p + geom_text(
    data = label_df,
    aes(x = x, y = y, label = block),
    inherit.aes = FALSE,
    hjust = 0,
    fontface = "bold",
    size = 3.2
  )
  
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  print(p)
  dev.off()
  
  invisible(list(
    plot = p,
    top_markers = top_markers,
    gene_layout = gene_layout
  ))
}

# =========================================================
# 4. Define hierarchical orders
# =========================================================

# Lymphoid: B/plasma -> CD4 states -> Treg -> CD8 states -> innate lymphoid/NK/NKT
lymphoid_order <- c(
  "B", "plasma",
  "CD4+ Tcm", "CD4+ Tex", "Treg",
  "CD8+ Tcm", "CD8 Tem", "CD8+ Teff", "CD8+ Tex",
  "NK", "NKT"
)

# Vasculature: arterial -> capillary/glomerular -> venous/mitotic
vasculature_order <- c(
  "pre-artery endothelial cells",
  "arteriole",
  "glomerulus",
  "lymphatic capillary",
  "venule",
  "mitotic venule"
)

# Myeloid: DC lineage -> pDC/tDC -> monocyte -> macrophage continuum -> mast
myeloid_order <- c(
  "cDC1", "cDC2", "pDC", "tDC",
  "monocyte",
  "M2a-like", "M2c-like", "M2d-like",
  "mast"
)

# =========================================================
# 5. Make the three PDFs
# =========================================================

lymphoid_markers <- compute_and_save_markers(
  obj = lymphoid,
  group_col = "cell_type",
  out_rds = file.path(outdir, "lymphoid_markers.rds"),
  out_csv = file.path(outdir, "lymphoid_markers.csv")
)

vasculature_markers <- compute_and_save_markers(
  obj = vasculature,
  group_col = "cell_type",
  out_rds = file.path(outdir, "vasculature_markers.rds"),
  out_csv = file.path(outdir, "vasculature_markers.csv")
)

myeloid_markers <- compute_and_save_markers(
  obj = myeloid,
  group_col = "cell_type",
  out_rds = file.path(outdir, "myeloid_markers.rds"),
  out_csv = file.path(outdir, "myeloid_markers.csv")
)

lymphoid_markers <- readRDS(file.path(outdir, "lymphoid_markers.rds"))
vasculature_markers <- readRDS(file.path(outdir, "vasculature_markers.rds"))
myeloid_markers <- readRDS(file.path(outdir, "myeloid_markers.rds"))

plot_dotplot_from_markers(
  obj = lymphoid,
  markers = lymphoid_markers,
  group_col = "cell_type",
  group_order = lymphoid_order,
  top_n = 5,
  title = "Wilms lymphoid: top 5 markers per subgroup",
  pdf_file = file.path(outdir, "Wilms_lymphoid_top5_marker_dotplot_clean.pdf"),
  width = 9,
  height = 14
)

plot_dotplot_from_markers(
  obj = vasculature,
  markers = vasculature_markers,
  group_col = "cell_type",
  group_order = vasculature_order,
  top_n = 5,
  title = "Wilms vasculature: top 5 markers per subgroup",
  pdf_file = file.path(outdir, "Wilms_vasculature_top5_marker_dotplot_clean.pdf"),
  width = 8,
  height = 10
)

plot_dotplot_from_markers(
  obj = myeloid,
  markers = myeloid_markers,
  group_col = "cell_type",
  group_order = myeloid_order,
  top_n = 5,
  title = "Wilms myeloid: top 5 markers per subgroup",
  pdf_file = file.path(outdir, "Wilms_myeloid_top5_marker_dotplot_clean.pdf"),
  width = 9,
  height = 12
)

message("Done. Files written to: ", outdir)
