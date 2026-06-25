#!/usr/bin/env Rscript

# =============================================================================
# Script: 03_plot_malignant_umaps_composition_markers.R
# Project: Wilms tumor plasticity manuscript
# Figure: Figure 1 - malignant compartment UMAPs, composition, and markers
#
# Purpose:
#   Generate malignant-cell UMAPs, subcompartment composition summaries,
#   marker dot plots, and per-sample abundance statistics for Wilms tumor
#   malignant states: blastema, stroma, and epithelium.
#
# Input:
#   data/processed/integrated_compartment_cancer_processed.RData
#       Seurat object named `s_objs` with UMAP coordinates and metadata columns:
#       subcompartment, condition, orig.ident.
#
# Outputs:
#   figures/figure_01/malignant_umap_subcompartment_no_label.pdf
#   figures/figure_01/malignant_umap_subcompartment_label.pdf
#   figures/figure_01/malignant_umap_condition_no_label.pdf
#   figures/figure_01/malignant_umap_condition_label.pdf
#   figures/figure_01/malignant_stacked_bar_condition_by_subcompartment.pdf
#   figures/figure_01/malignant_stacked_bar_condition_by_subcompartment_labeled.pdf
#   figures/figure_01/malignant_dotplot_top10_markers_per_subcompartment.pdf
#   figures/figure_01/malignant_boxplot_subcompartment_fraction_by_condition.pdf
#
#   results/figure_01/all_markers_subcompartment.csv
#   results/figure_01/top10_markers_per_subcompartment.csv
#   results/figure_01/subcompartment_fraction_by_condition_stats.csv
#
# Run from repository root:
#   Rscript scripts/figure_01/03_plot_malignant_umaps_composition_markers.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(rstatix)
  library(scales)
})

# -----------------------------------------------------------------------------
# User-configurable paths
# -----------------------------------------------------------------------------
cancer_path <- file.path("data", "processed", "integrated_compartment_cancer_processed.RData")
out_dir <- file.path("figures", "figure_01")
results_dir <- file.path("results", "figure_01")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cancer_path)) {
  stop("Missing required input file: ", cancer_path, call. = FALSE)
}

load(cancer_path)
if (!exists("s_objs")) {
  stop("The input RData file must contain a Seurat object named `s_objs`.", call. = FALSE)
}

cancer <- s_objs

required_cols <- c("subcompartment", "condition", "orig.ident")
missing_cols <- setdiff(required_cols, colnames(cancer@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing columns in cancer@meta.data: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# Factor ordering and color palettes
# -----------------------------------------------------------------------------
cancer$subcompartment <- factor(cancer$subcompartment, levels = c("blastema", "stroma", "epithelium"))
cancer$condition <- factor(cancer$condition, levels = c("untreated", "neoadjuvant"))

subcompartment_cols <- c(
  "blastema" = "#E64B35",
  "stroma" = "#4DBBD5",
  "epithelium" = "#00A087"
)

condition_cols <- c(
  "untreated" = "#2166AC",
  "neoadjuvant" = "#B2182B"
)

# -----------------------------------------------------------------------------
# Helper functions
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

save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  out_file <- file.path(out_dir, filename)
  ggsave(out_file, plot_obj, width = width, height = height)
  message("Saved: ", out_file)
}

# -----------------------------------------------------------------------------
# 1. Malignant UMAP panels
# -----------------------------------------------------------------------------
p_sub_no_label <- plot_umap(cancer, "subcompartment", subcompartment_cols, "UMAP - Malignant subcompartment", label = FALSE)
p_sub_label <- plot_umap(cancer, "subcompartment", subcompartment_cols, "UMAP - Malignant subcompartment", label = TRUE)
p_cond_no_label <- plot_umap(cancer, "condition", condition_cols, "UMAP - Condition", label = FALSE)
p_cond_label <- plot_umap(cancer, "condition", condition_cols, "UMAP - Condition", label = TRUE)

save_plot(p_sub_no_label, "malignant_umap_subcompartment_no_label.pdf")
save_plot(p_sub_label, "malignant_umap_subcompartment_label.pdf")
save_plot(p_cond_no_label, "malignant_umap_condition_no_label.pdf")
save_plot(p_cond_label, "malignant_umap_condition_label.pdf")

# -----------------------------------------------------------------------------
# 2. Condition composition within each malignant subcompartment
# -----------------------------------------------------------------------------
bar_df <- cancer@meta.data %>%
  dplyr::count(subcompartment, condition, name = "n") %>%
  dplyr::group_by(subcompartment) %>%
  dplyr::mutate(
    percent = 100 * n / sum(n),
    percent_label = paste0(round(percent, 1), "%")
  ) %>%
  dplyr::ungroup()

p_bar_no_labels <- ggplot(bar_df, aes(x = subcompartment, y = percent, fill = condition)) +
  geom_bar(stat = "identity", width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "Condition composition within each malignant subcompartment",
    x = "Malignant subcompartment",
    y = "Percentage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_bar_labels <- p_bar_no_labels +
  geom_text(
    aes(label = percent_label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white",
    fontface = "bold"
  )

save_plot(p_bar_no_labels, "malignant_stacked_bar_condition_by_subcompartment.pdf", width = 6, height = 5)
save_plot(p_bar_labels, "malignant_stacked_bar_condition_by_subcompartment_labeled.pdf", width = 6, height = 5)

# -----------------------------------------------------------------------------
# 3. Top marker dot plot per malignant subcompartment
# -----------------------------------------------------------------------------
Idents(cancer) <- cancer$subcompartment
DefaultAssay(cancer) <- "RNA"

markers_sub <- FindAllMarkers(
  cancer,
  only.pos = TRUE,
  assay = "RNA",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10_sub <- markers_sub %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

write.csv(markers_sub, file.path(results_dir, "all_markers_subcompartment.csv"), row.names = FALSE)
write.csv(top10_sub, file.path(results_dir, "top10_markers_per_subcompartment.csv"), row.names = FALSE)

avg_exp <- AverageExpression(
  cancer,
  assays = "RNA",
  group.by = "subcompartment",
  slot = "data",
  return.seurat = FALSE
)

avg_mat <- avg_exp$RNA
sel_genes <- unique(top10_sub$gene)
avg_mat2 <- avg_mat[intersect(sel_genes, rownames(avg_mat)), , drop = FALSE]

# Hierarchical ordering of malignant subcompartments by average marker expression.
hc_groups <- hclust(dist(t(avg_mat2)))
group_order <- hc_groups$labels[hc_groups$order]

cancer$subcompartment <- factor(cancer$subcompartment, levels = group_order)
Idents(cancer) <- cancer$subcompartment

top10_sub$cluster <- factor(top10_sub$cluster, levels = group_order)
gene_order_df <- top10_sub %>%
  dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  dplyr::distinct(gene, cluster, avg_log2FC)

gene_order <- rev(gene_order_df$gene)

p_dot <- DotPlot(
  cancer,
  features = gene_order,
  group.by = "subcompartment",
  assay = "RNA"
)

# Force dot outlines to black while preserving expression-based fill scale.
p_dot$layers[[1]]$aes_params$colour <- "black"

p_dot <- p_dot +
  coord_flip() +
  labs(
    title = "Top 10 markers per malignant subcompartment",
    x = "Malignant subcompartment",
    y = "Markers"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

save_plot(p_dot, "malignant_dotplot_top10_markers_per_subcompartment.pdf", width = 5, height = 8)

# -----------------------------------------------------------------------------
# 4. Per-sample malignant subcompartment abundance by condition
# -----------------------------------------------------------------------------
# Restore the biologically intuitive order for abundance plotting.
cancer$subcompartment <- factor(cancer$subcompartment, levels = c("blastema", "stroma", "epithelium"))

frac_df <- cancer@meta.data %>%
  dplyr::select(orig.ident, condition, subcompartment) %>%
  dplyr::count(orig.ident, condition, subcompartment, name = "n_cells") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(
    total_cells_sample = sum(n_cells),
    fraction = n_cells / total_cells_sample
  ) %>%
  dplyr::ungroup()

stat_df <- frac_df %>%
  dplyr::group_by(subcompartment) %>%
  rstatix::wilcox_test(fraction ~ condition) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::ungroup()

ypos_df <- frac_df %>%
  dplyr::group_by(subcompartment) %>%
  dplyr::summarise(y.position = max(fraction, na.rm = TRUE) * 1.12, .groups = "drop")

stat_df <- stat_df %>%
  dplyr::left_join(ypos_df, by = "subcompartment") %>%
  dplyr::mutate(
    p_label = paste0("BH p = ", signif(p.adj, 3)),
    xnum = match(subcompartment, levels(frac_df$subcompartment)),
    xmin = xnum - 0.2,
    xmax = xnum + 0.2,
    xmid = xnum
  )

write.csv(stat_df, file.path(results_dir, "subcompartment_fraction_by_condition_stats.csv"), row.names = FALSE)

p_box <- ggplot(frac_df, aes(x = subcompartment, y = fraction, fill = condition)) +
  geom_boxplot(
    aes(color = condition),
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.35,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(color = condition),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    size = 2,
    alpha = 0.85
  ) +
  geom_segment(
    data = stat_df,
    aes(x = xmin, xend = xmax, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = stat_df,
    aes(x = xmin, xend = xmin, y = y.position * 0.985, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_segment(
    data = stat_df,
    aes(x = xmax, xend = xmax, y = y.position * 0.985, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.4
  ) +
  geom_text(
    data = stat_df,
    aes(x = xmid, y = y.position * 1.02, label = p_label),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = condition_cols) +
  scale_color_manual(values = condition_cols) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.18))
  ) +
  labs(
    title = "Malignant subcompartment abundance by condition",
    x = "Malignant subcompartment",
    y = "Fraction of malignant cells per sample"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

save_plot(p_box, "malignant_boxplot_subcompartment_fraction_by_condition.pdf", width = 8, height = 6)

message("Finished malignant UMAP, composition, marker, and abundance analyses.")
