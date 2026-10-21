# ============================================================
# Figure 2 malignant-cell UMAPs, subcompartment composition, marker dot plots, and abundance statistics
# ============================================================
#
# GitHub script name : 01_plot_malignant_umaps_composition_markers.R
# Original file      : Fig.2.wilms-script.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Loads malignant/cancer Seurat object and generates Figure 2-level malignant compartment summaries.
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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(rstatix)
  library(scales)
})

# =========================================================
# 0. Setup
# =========================================================

outdir <- "~/Documents/Izar_Group/wilms/Figures/v2/Fig.2.v2"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

obj_path <- "~/Documents/Izar_Group/wilms/integrated_compartment_cancer_processed.RData"

load(obj_path)   # loads s_objs
cancer <- s_objs

DefaultAssay(cancer) <- "RNA"

# Required metadata columns
required_cols <- c("subcompartment", "condition", "orig.ident")
missing_cols <- setdiff(required_cols, colnames(cancer@meta.data))

if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# Factor order
subcompartment_order <- c("blastema", "stroma", "epithelium")
condition_order <- c("untreated", "neoadjuvant")

cancer$subcompartment <- factor(cancer$subcompartment, levels = subcompartment_order)
cancer$condition <- factor(cancer$condition, levels = condition_order)

# Colors
subcompartment_cols <- c(
  "blastema"   = "#E64B35",
  "stroma"     = "#4DBBD5",
  "epithelium" = "#00A087"
)

condition_cols <- c(
  "untreated"   = "#7F7F7F",
  "neoadjuvant" = "#E69F00"
)

# =========================================================
# 1. UMAP plots
# =========================================================

make_umap <- function(obj, group_col, cols, label = FALSE, title = NULL) {
  DimPlot(
    obj,
    reduction = "umap",
    group.by = group_col,
    label = label,
    repel = TRUE,
    label.size = 5,
    pt.size = 0.1,
    raster = FALSE
  ) +
    scale_color_manual(values = cols, na.value = "grey80") +
    ggtitle(title) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank()
    )
}

p_sub_no_label <- make_umap(
  cancer,
  group_col = "subcompartment",
  cols = subcompartment_cols,
  label = FALSE,
  title = "UMAP - Subcompartment"
)

p_sub_label <- make_umap(
  cancer,
  group_col = "subcompartment",
  cols = subcompartment_cols,
  label = TRUE,
  title = "UMAP - Subcompartment"
)

p_cond_no_label <- make_umap(
  cancer,
  group_col = "condition",
  cols = condition_cols,
  label = FALSE,
  title = "UMAP - Condition"
)

p_cond_label <- make_umap(
  cancer,
  group_col = "condition",
  cols = condition_cols,
  label = TRUE,
  title = "UMAP - Condition"
)

ggsave(file.path(outdir, "UMAP_subcompartment_no_label.pdf"), p_sub_no_label, width = 8, height = 6)
ggsave(file.path(outdir, "UMAP_subcompartment_label.pdf"),    p_sub_label,    width = 8, height = 6)
ggsave(file.path(outdir, "UMAP_condition_no_label.pdf"),      p_cond_no_label, width = 8, height = 6)
ggsave(file.path(outdir, "UMAP_condition_label.pdf"),         p_cond_label,    width = 8, height = 6)

# =========================================================
# 2. Stacked bar plots: condition within subcompartment
# =========================================================

bar_df <- cancer@meta.data %>%
  count(subcompartment, condition, name = "n") %>%
  group_by(subcompartment) %>%
  mutate(
    percent = 100 * n / sum(n),
    percent_label = paste0(round(percent, 1), "%")
  ) %>%
  ungroup()

p_bar_percent <- ggplot(bar_df, aes(x = subcompartment, y = percent, fill = condition)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "Condition composition within each subcompartment",
    x = "Subcompartment",
    y = "Percentage"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_bar_percent_labels <- p_bar_percent +
  geom_text(
    aes(label = percent_label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white",
    fontface = "bold"
  )

p_bar_counts <- ggplot(bar_df, aes(x = subcompartment, y = n, fill = condition)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "Cell counts by condition within each subcompartment",
    x = "Subcompartment",
    y = "Cell count"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

ggsave(file.path(outdir, "stacked_bar_condition_percent_no_labels.pdf"), p_bar_percent, width = 6, height = 5)
ggsave(file.path(outdir, "stacked_bar_condition_percent_with_labels.pdf"), p_bar_percent_labels, width = 6, height = 5)
ggsave(file.path(outdir, "stacked_bar_condition_counts.pdf"), p_bar_counts, width = 6, height = 5)

# =========================================================
# 3. Marker discovery and black-dot marker DotPlot
# =========================================================

Idents(cancer) <- cancer$subcompartment

markers_sub <- FindAllMarkers(
  cancer,
  only.pos = TRUE,
  assay = "RNA",
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10_sub <- markers_sub %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(markers_sub, file.path(outdir, "all_markers_subcompartment.csv"), row.names = FALSE)
write.csv(top10_sub, file.path(outdir, "top10_markers_per_subcompartment.csv"), row.names = FALSE)

sel_genes <- unique(top10_sub$gene)

avg_exp <- AverageExpression(
  cancer,
  assays = "RNA",
  group.by = "subcompartment",
  layer = "data",
  return.seurat = FALSE
)

avg_mat <- avg_exp$RNA
avg_mat2 <- avg_mat[intersect(sel_genes, rownames(avg_mat)), , drop = FALSE]

if (nrow(avg_mat2) >= 2 && ncol(avg_mat2) >= 2) {
  hc_groups <- hclust(dist(t(avg_mat2)))
  group_order <- hc_groups$labels[hc_groups$order]
} else {
  group_order <- subcompartment_order
}

cancer$subcompartment <- factor(cancer$subcompartment, levels = group_order)
Idents(cancer) <- cancer$subcompartment

top10_sub$cluster <- factor(top10_sub$cluster, levels = group_order)

gene_order <- top10_sub %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  distinct(gene, cluster, avg_log2FC) %>%
  pull(gene) %>%
  rev()

p_dot <- DotPlot(
  cancer,
  features = gene_order,
  group.by = "subcompartment",
  assay = "RNA"
)

p_dot$layers[[1]]$aes_params$colour <- "black"
  
p_dot <- p_dot +
  coord_flip() +
  labs(
    title = "Top 10 markers per subcompartment",
    x = "Subcompartment",
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

ggsave(file.path(outdir, "dotplot_top10_markers_per_subcompartment_black.pdf"), p_dot, width = 5, height = 8)


# =========================================================
# 4. Per-sample subcompartment abundance by condition
# =========================================================

frac_df <- cancer@meta.data %>%
  select(orig.ident, condition, subcompartment) %>%
  count(orig.ident, condition, subcompartment, name = "n_cells") %>%
  group_by(orig.ident) %>%
  mutate(
    total_cells_sample = sum(n_cells),
    fraction = n_cells / total_cells_sample
  ) %>%
  ungroup()

stat_df <- frac_df %>%
  group_by(subcompartment) %>%
  wilcox_test(fraction ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

ypos_df <- frac_df %>%
  group_by(subcompartment) %>%
  summarise(y.position = max(fraction, na.rm = TRUE) * 1.12, .groups = "drop")

stat_df <- stat_df %>%
  left_join(ypos_df, by = "subcompartment") %>%
  mutate(
    p_label = paste0(p.adj.signif, " (BH p = ", signif(p.adj, 3), ")"),
    xnum = match(subcompartment, levels(frac_df$subcompartment)),
    xmin = xnum - 0.2,
    xmax = xnum + 0.2,
    xmid = xnum
  )

write.csv(frac_df, file.path(outdir, "sample_fraction_subcompartment_by_condition.csv"), row.names = FALSE)
write.csv(stat_df, file.path(outdir, "sample_fraction_subcompartment_stats.csv"), row.names = FALSE)

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
    aes(x = xmid, y = y.position * 1.03, label = p_label),
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
    title = "Subcompartment abundance by condition",
    x = "Subcompartment",
    y = "Fraction of cells per sample"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_blank()
  )

ggsave(file.path(outdir, "boxplot_subcompartment_fraction_by_condition_with_pvalues.pdf"), p_box, width = 8, height = 6)

# =========================================================
# Done
# =========================================================

message("Done. All outputs saved to: ", normalizePath(outdir))
