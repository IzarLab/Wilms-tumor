# ============================================================
# Figure 5 TME lymphoid, myeloid, and vasculature composition
# ============================================================
#
# GitHub script name : 01_plot_tme_lymphoid_myeloid_vasculature_composition.R
# Original file      : Fig.5.wilms-TME.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Generates TME composition panels for lymphoid, myeloid, and vascular compartments with aligned clinical annotations.
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

# ============================================================
# Wilms TME: lymphoid / myeloid / vasculature
# Correct 100% stacked plot + aligned clinical table
# W004 removed before plotting because no selected cells
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# -----------------------------
# Input / output
# -----------------------------
infile <- "~/Documents/Izar_Group/wilms/TME.RData"
outdir <- "~/Documents/Izar_Group/wilms/TME_lymphoid_myeloid_vasculature"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

load(infile)
obj <- tme_

# ============================================================
# Sample metadata
# ============================================================

sample_meta <- tibble::tribble(
  ~sample, ~tumor_analyzed_code, ~neoadjuvant_code, ~response_code, ~tumor_histology_code, ~disease_stage, ~disease_recurrence_code,
  "W612", 1, "N/A", "N/A", 1, 2, "1",
  "W340", 1, "N/A", "N/A", 1, 4, "1",
  "W011", 1, "N/A", "N/A", 1, 3, "2",
  "W132", 1, "N/A", "N/A", 1, 4, "2",
  "W495", 1, "N/A", "N/A", 1, 2, "2",
  "W008", 1, "N/A", "N/A", 1, 3, "2",
  "W007", 1, "N/A", "N/A", 1, 2, "2",
  "W005", 1, "N/A", "N/A", 1, 2, "2",
  "W004", 1, "N/A", "N/A", 1, 1, "1",
  "W241", 1, "N/A", "N/A", 2, 3, "UNK",
  "W040", 1, "N/A", "N/A", 1, 2, "1",
  "W012", 1, "N/A", "N/A", 1, 3, "2",
  "W002", 1, "3",   "1",   2, 4, "1",
  "W050", 1, "5*",  "2",   1, 3, "1",
  "W010", 1, "1",   "1",   1, 2, "2",
  "W944", 1, "2",   "1",   1, 4, "2",
  "W009", 1, "1",   "3",   1, 5, "2",
  "W184", 1, "2",   "1",   2, 3, "2",
  "W180", 2, "2",   "2",   1, 3, "1",
  "W370", 1, "2",   "1",   2, 5, "2",
  "W006", 2, "2",   "UNK", 1, 3, "1",
  "W003", 1, "1",   "1",   1, 4, "1",
  "W441", 1, "1",   "2",   1, 4, "2",
  "W001", 2, "UNK", "UNK", 1, 1, "1"
)

sample_meta <- sample_meta %>%
  mutate(
    tumor_analyzed = case_when(
      tumor_analyzed_code == 1 ~ "Primary",
      tumor_analyzed_code == 2 ~ "Recurrence",
      TRUE ~ as.character(tumor_analyzed_code)
    ),
    neoadjuvant = case_when(
      neoadjuvant_code == "1" ~ "EE-4A",
      neoadjuvant_code == "2" ~ "DD-4A",
      neoadjuvant_code == "3" ~ "I",
      neoadjuvant_code == "4" ~ "M",
      neoadjuvant_code %in% c("5", "5*") ~ "V+D",
      neoadjuvant_code == "N/A" ~ "Untreated",
      neoadjuvant_code == "UNK" ~ "UNK",
      TRUE ~ neoadjuvant_code
    ),
    response = case_when(
      response_code == "1" ~ "Dec",
      response_code == "2" ~ "Stable",
      response_code == "3" ~ "Inc",
      response_code == "N/A" ~ "N/A",
      response_code == "UNK" ~ "UNK",
      TRUE ~ response_code
    ),
    histology = case_when(
      tumor_histology_code == 1 ~ "Favorable",
      tumor_histology_code == 2 ~ "Anaplastic",
      TRUE ~ as.character(tumor_histology_code)
    ),
    disease_stage = as.character(disease_stage),
    recurrence = case_when(
      disease_recurrence_code == "1" ~ "Yes",
      disease_recurrence_code == "2" ~ "No",
      disease_recurrence_code == "UNK" ~ "UNK",
      TRUE ~ disease_recurrence_code
    )
  )

# ============================================================
# Subset selected compartments
# ============================================================

keep_subcomp <- c("vasculature", "myeloid", "lymphoid")

subcomp_cols <- c(
  "vasculature" = "#1B9E77",
  "myeloid"     = "#D95F02",
  "lymphoid"    = "#7570B3"
)

obj.sub <- subset(
  obj,
  subset = subcompartment %in% keep_subcomp
)

obj.sub$subcompartment <- factor(
  obj.sub$subcompartment,
  levels = keep_subcomp
)

obj.sub$sample_plot <- as.character(obj.sub$orig.ident)

obj.sub$sample_plot <- case_when(
  obj.sub$sample_plot == "W1"  ~ "W001",
  obj.sub$sample_plot == "W2"  ~ "W002",
  obj.sub$sample_plot == "W3"  ~ "W003",
  obj.sub$sample_plot == "W4"  ~ "W004",
  obj.sub$sample_plot == "W5"  ~ "W005",
  obj.sub$sample_plot == "W6"  ~ "W006",
  obj.sub$sample_plot == "W7"  ~ "W007",
  obj.sub$sample_plot == "W8"  ~ "W008",
  obj.sub$sample_plot == "W9"  ~ "W009",
  obj.sub$sample_plot == "W10" ~ "W010",
  obj.sub$sample_plot == "W11" ~ "W011",
  obj.sub$sample_plot == "W12" ~ "W012",
  obj.sub$sample_plot == "W40" ~ "W040",
  obj.sub$sample_plot == "W50" ~ "W050",
  TRUE ~ obj.sub$sample_plot
)

# ============================================================
# Re-normalize from raw RNA counts
# ============================================================

DefaultAssay(obj.sub) <- "RNA"

obj.sub <- NormalizeData(obj.sub)
obj.sub <- FindVariableFeatures(obj.sub, selection.method = "vst", nfeatures = 3000)
obj.sub <- ScaleData(obj.sub, features = VariableFeatures(obj.sub))
obj.sub <- RunPCA(obj.sub, features = VariableFeatures(obj.sub), npcs = 50)
obj.sub <- FindNeighbors(obj.sub, dims = 1:30)
obj.sub <- FindClusters(obj.sub, resolution = 0.4)
obj.sub <- RunUMAP(
  obj.sub,
  dims = 1:30,
  reduction.name = "umap.recomputed",
  reduction.key = "UMAPre_"
)

saveRDS(
  obj.sub,
  file.path(outdir, "TME_Lymphoid_Myeloid_Vasculature_processed.rds")
)

# ============================================================
# UMAP
# ============================================================

p_umap <- DimPlot(
  obj.sub,
  reduction = "umap.recomputed",
  group.by = "subcompartment",
  cols = subcomp_cols,
  pt.size = 0.15,
  raster = FALSE
) +
  theme_classic(base_size = 14)

ggsave(
  file.path(outdir, "UMAP_lymphoid_myeloid_vasculature.pdf"),
  p_umap,
  width = 7,
  height = 6,
  device = cairo_pdf
)

# ============================================================
# Correct 100% stack data
# Remove W004 before plotting
# ============================================================

sample_order_full <- sample_meta$sample

count_df <- obj.sub@meta.data %>%
  filter(
    !is.na(sample_plot),
    sample_plot %in% sample_order_full,
    subcompartment %in% keep_subcomp
  ) %>%
  count(sample_plot, subcompartment, name = "n")

total_df <- count_df %>%
  group_by(sample_plot) %>%
  summarise(total_selected_cells = sum(n), .groups = "drop")

sample_order_plot <- sample_order_full[
  sample_order_full %in% total_df$sample_plot
]

# Explicitly remove W004 because it has no cells in this subset
sample_order_plot <- setdiff(sample_order_plot, "W004")

prop_df <- count_df %>%
  filter(sample_plot %in% sample_order_plot) %>%
  complete(
    sample_plot = factor(sample_order_plot, levels = sample_order_plot),
    subcompartment = factor(keep_subcomp, levels = keep_subcomp),
    fill = list(n = 0)
  ) %>%
  left_join(total_df, by = "sample_plot") %>%
  mutate(
    prop = 100 * n / total_selected_cells
  ) %>%
  left_join(sample_meta, by = c("sample_plot" = "sample"))

prop_df$sample_plot <- factor(prop_df$sample_plot, levels = sample_order_plot)
prop_df$subcompartment <- factor(prop_df$subcompartment, levels = keep_subcomp)

# Manual ymin/ymax prevents blank stack segments
plot_df <- prop_df %>%
  arrange(sample_plot, subcompartment) %>%
  group_by(sample_plot) %>%
  mutate(
    ymin = cumsum(lag(prop, default = 0)),
    ymax = ymin + prop,
    ymid = (ymin + ymax) / 2
  ) %>%
  ungroup()

write.csv(
  plot_df,
  file.path(outdir, "stackplot_values_100percent_no_W004.csv"),
  row.names = FALSE
)

# ============================================================
# Aligned table below plot using same x-axis positions
# ============================================================

table_df <- sample_meta %>%
  filter(sample %in% sample_order_plot) %>%
  mutate(sample = factor(sample, levels = sample_order_plot)) %>%
  arrange(sample) %>%
  select(
    sample,
    tumor_analyzed,
    neoadjuvant,
    response,
    histology,
    disease_stage,
    recurrence
  ) %>%
  rename(
    Sample = sample,
    Tumor = tumor_analyzed,
    Treatment = neoadjuvant,
    Response = response,
    Histology = histology,
    Stage = disease_stage,
    Recurrence = recurrence
  ) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "row",
    values_to = "value"
  )

table_rows <- c(
  "Sample",
  "Tumor",
  "Treatment",
  "Response",
  "Histology",
  "Stage",
  "Recurrence"
)

table_y <- tibble(
  row = table_rows,
  y = c(-7, -15, -23, -31, -39, -47, -55)
)

table_plot_df <- table_df %>%
  bind_rows(
    tibble(
      Sample = factor(sample_order_plot, levels = sample_order_plot),
      row = "Sample",
      value = sample_order_plot
    )
  ) %>%
  left_join(table_y, by = "row") %>%
  mutate(
    row = factor(row, levels = table_rows),
    Sample = factor(Sample, levels = sample_order_plot)
  )

row_label_df <- table_y %>%
  mutate(
    x = 0.25,
    label = row
  )

# ============================================================
# Final stack plot with aligned table
# ============================================================

# p_stack_table <- ggplot() +
#   geom_rect(
#     data = plot_df,
#     aes(
#       xmin = as.numeric(sample_plot) - 0.42,
#       xmax = as.numeric(sample_plot) + 0.42,
#       ymin = ymin,
#       ymax = ymax,
#       fill = subcompartment
#     ),
#     color = "black",
#     linewidth = 0.15
#   ) +
#   geom_text(
#     data = plot_df %>% filter(prop >= 4),
#     aes(
#       x = as.numeric(sample_plot),
#       y = ymid,
#       label = round(prop, 1)
#     ),
#     size = 3,
#     fontface = "bold",
#     color = "black"
#   ) +
#   geom_text(
#     data = table_plot_df,
#     aes(
#       x = as.numeric(Sample),
#       y = y,
#       label = value
#     ),
#     size = 2.4,
#     fontface = "bold",
#     color = "black"
#   ) +
#   geom_text(
#     data = row_label_df,
#     aes(
#       x = x,
#       y = y,
#       label = label
#     ),
#     size = 2.8,
#     fontface = "bold",
#     hjust = 1,
#     color = "black"
#   ) +
#   geom_hline(yintercept = -2.5, linewidth = 0.3) +
#   scale_fill_manual(
#     values = subcomp_cols,
#     breaks = c("lymphoid", "myeloid", "vasculature"),
#     name = "Compartment"
#   ) +
#   scale_x_continuous(
#     breaks = seq_along(sample_order_plot),
#     labels = sample_order_plot,
#     limits = c(0.3, length(sample_order_plot) + 0.8),
#     expand = c(0, 0)
#   ) +
#   scale_y_continuous(
#     breaks = c(0, 25, 50, 75, 100),
#     limits = c(-60, 102),
#     expand = c(0, 0)
#   ) +
#   labs(
#     x = NULL,
#     y = "Proportion within lymphoid + myeloid + vasculature (%)"
#   ) +
#   coord_cartesian(clip = "off") +
#   theme_classic(base_size = 13) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.line.x = element_blank(),
#     axis.title.y = element_text(face = "bold"),
#     axis.text.y = element_text(face = "bold"),
#     legend.title = element_text(face = "bold"),
#     legend.text = element_text(face = "bold"),
#     plot.margin = margin(10, 20, 10, 50)
#   )
# 
# ggsave(
#   file.path(outdir, "StackedBar_100percent_aligned_table_no_W004.pdf"),
#   p_stack_table,
#   width = 16,
#   height = 9,
#   device = "pdf"
# )
# 

# ============================================================
# Table geometry
# ============================================================

table_rows <- c(
  "Sample",
  "Tumor",
  "Treatment",
  "Response",
  "Histology",
  "Stage",
  "Recurrence"
)

table_y <- tibble(
  row = table_rows,
  y = c(-7, -15, -23, -31, -39, -47, -55)
)

# -----------------------------
# Table values
# -----------------------------
table_df <- sample_meta %>%
  filter(sample %in% sample_order_plot) %>%
  mutate(sample = factor(sample, levels = sample_order_plot)) %>%
  arrange(sample) %>%
  select(
    sample,
    tumor_analyzed,
    neoadjuvant,
    response,
    histology,
    disease_stage,
    recurrence
  ) %>%
  rename(
    Sample = sample,
    Tumor = tumor_analyzed,
    Treatment = neoadjuvant,
    Response = response,
    Histology = histology,
    Stage = disease_stage,
    Recurrence = recurrence
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "row",
    values_to = "value"
  )

table_df <- table_df %>%
  left_join(table_y, by = "row")

# -----------------------------
# Numeric x positions
# -----------------------------
sample_x_df <- tibble(
  sample = sample_order_plot,
  xpos = seq_along(sample_order_plot)
)

table_df <- table_df %>%
  rename(sample = value) %>%
  mutate(
    sample_original = sample
  )

# rebuild table correctly
table_df <- sample_meta %>%
  filter(sample %in% sample_order_plot) %>%
  mutate(sample = factor(sample, levels = sample_order_plot)) %>%
  arrange(sample)

table_long <- bind_rows(
  
  tibble(
    sample = sample_order_plot,
    row = "Sample",
    value = sample_order_plot
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Tumor",
    value = table_df$tumor_analyzed
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Treatment",
    value = table_df$neoadjuvant
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Response",
    value = table_df$response
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Histology",
    value = table_df$histology
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Stage",
    value = table_df$disease_stage
  ),
  
  tibble(
    sample = table_df$sample,
    row = "Recurrence",
    value = table_df$recurrence
  )
  
)

table_long <- table_long %>%
  left_join(table_y, by = "row") %>%
  left_join(
    sample_x_df,
    by = c("sample" = "sample")
  )

# ============================================================
# Borders
# ============================================================

# Vertical borders
vertical_lines <- tibble(
  x = seq(0.5, length(sample_order_plot) + 0.5, by = 1)
)

# Horizontal borders
horizontal_lines <- tibble(
  y = c(-3, -11, -19, -27, -35, -43, -51, -59)
)

# Row labels
row_label_df <- tibble(
  row = table_rows,
  y = c(-7, -15, -23, -31, -39, -47, -55),
  label = c(
    "Sample",
    "Tumor",
    "Treatment",
    "Response",
    "Histology",
    "Stage",
    "Recurrence"
  )
)

# ============================================================
# Final combined plot
# ============================================================

p_stack_table <- ggplot() +
  
  # -----------------------------
# STACKED BARS
# -----------------------------
geom_rect(
  data = plot_df,
  aes(
    xmin = as.numeric(sample_plot) - 0.42,
    xmax = as.numeric(sample_plot) + 0.42,
    ymin = ymin,
    ymax = ymax,
    fill = subcompartment
  ),
  color = "black",
  linewidth = 0.15
) +
  
  geom_text(
    data = plot_df %>% filter(prop >= 4),
    aes(
      x = as.numeric(sample_plot),
      y = ymid,
      label = round(prop, 1)
    ),
    size = 3,
    fontface = "bold"
  ) +
  
  # -----------------------------
# TABLE BORDERS
# -----------------------------
geom_vline(
  data = vertical_lines,
  aes(xintercept = x),
  ymin = 0,
  ymax = 1,
  linewidth = 0.25,
  color = "grey50"
) +
  
  geom_hline(
    data = horizontal_lines,
    aes(yintercept = y),
    linewidth = 0.25,
    color = "grey50"
  ) +
  
  # -----------------------------
# TABLE TEXT
# -----------------------------
geom_text(
  data = table_long,
  aes(
    x = xpos,
    y = y,
    label = value
  ),
  size = 2.5,
  fontface = "bold"
) +
  
  # -----------------------------
# ROW LABELS
# -----------------------------
geom_text(
  data = row_label_df,
  aes(
    x = 0.15,
    y = y,
    label = label
  ),
  hjust = 1,
  size = 2.8,
  fontface = "bold"
) +
  
  # -----------------------------
# TABLE LEFT BORDER
# -----------------------------
geom_vline(
  xintercept = 0.5,
  linewidth = 0.35,
  color = "black"
) +
  
  # -----------------------------
# TABLE TOP BORDER
# -----------------------------
geom_hline(
  yintercept = -3,
  linewidth = 0.35,
  color = "black"
) +
  
  scale_fill_manual(
    values = subcomp_cols,
    breaks = c("lymphoid", "myeloid", "vasculature"),
    name = "Compartment"
  ) +
  
  scale_x_continuous(
    breaks = seq_along(sample_order_plot),
    labels = sample_order_plot,
    limits = c(-0.8, length(sample_order_plot) + 0.6),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(-60, 102),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Proportion within lymphoid + myeloid + vasculature (%)"
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    
    axis.title.y = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    
    plot.margin = margin(10, 20, 10, 60)
  )

# ============================================================
# SAVE
# ============================================================

ggsave(
  file.path(outdir, "StackedBar_100percent_aligned_table_bordered.pdf"),
  p_stack_table,
  width = 16,
  height = 9,
  device = "pdf"
)


message("Done. Outputs written to: ", outdir)


