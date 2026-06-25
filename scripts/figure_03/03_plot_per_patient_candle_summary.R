# ============================================================
# Per-patient CANDLE plasticity summary
# ============================================================
#
# GitHub script name : 03_plot_per_patient_candle_summary.R
# Original file      : per-patient-CANDLE-figure.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Generates per-patient CANDLE switching/selection summary panels and includes optional sample-level CNV/Numbat helper code.
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
#   results/figure_03/
#
# ============================================================

# =========================================================
# Per-patient CANDLE summary panel for Wilms data
# Input: wilms_wmNMF_K25.plasticity.per_patient.tsv
# Output:
#   1) wilms_CANDLE_per_patient_summary_panel.pdf
#   2) wilms_CANDLE_per_patient_summary_panel.png
# =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(scales)
})

# ----------------------------
# 1. Read input
# ----------------------------
infile <- "wilms_wmNMF_K25.plasticity.per_patient.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

# Expected columns:
# patient, within_switching, between_selection, total_change,
# switch_fraction, condition, histology, n_cells

# ----------------------------
# 2. Basic checks
# ----------------------------
req_cols <- c(
  "patient",
  "within_switching",
  "between_selection",
  "total_change",
  "switch_fraction",
  "condition",
  "histology",
  "n_cells"
)

missing_cols <- setdiff(req_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Clean up
df <- df %>%
  mutate(
    patient = as.character(patient),
    condition = as.character(condition),
    histology = as.character(histology),
    n_cells = as.numeric(n_cells),
    within_switching = as.numeric(within_switching),
    between_selection = as.numeric(between_selection),
    total_change = as.numeric(total_change),
    switch_fraction = as.numeric(switch_fraction)
  )

# Optional derived quantity
df <- df %>%
  mutate(
    between_fraction = 1 - switch_fraction
  )

# ----------------------------
# 3. Order patients
#    You can change ordering to total_change if preferred
# ----------------------------
df_ord <- df %>%
  arrange(desc(within_switching), desc(total_change)) %>%
  mutate(patient = factor(patient, levels = patient))

# Long format for panel A
df_long <- df_ord %>%
  select(patient, within_switching, between_selection) %>%
  pivot_longer(
    cols = c(within_switching, between_selection),
    names_to = "component",
    values_to = "score"
  ) %>%
  mutate(
    component = recode(
      component,
      "within_switching" = "Within switching",
      "between_selection" = "Between selection"
    ),
    component = factor(component, levels = c("Within switching", "Between selection"))
  )

# ----------------------------
# 4. Theme
# ----------------------------
theme_candle <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 12),
    plot.title.position = "plot",
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    strip.background = element_rect(fill = "white", color = "black")
  )

# ----------------------------
# 5. Panel A:
#    Per-patient decomposition
# ----------------------------
pA <- ggplot(df_long, aes(x = patient, y = score, fill = component)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.68,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "Within switching" = "grey25",
      "Between selection" = "grey75"
    )
  ) +
  labs(
    title = "A. Per-patient CANDLE decomposition",
    x = NULL,
    y = "CANDLE score",
    fill = NULL
  ) +
  theme_candle +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.y = element_line(color = "grey90")
  )

# ----------------------------
# 6. Panel B:
#    Within vs between scatter
# ----------------------------
# size scaling for n_cells
# Adjust range if points are too large/small
pB <- ggplot(
  df,
  aes(
    x = within_switching,
    y = between_selection,
    color = histology,
    shape = condition,
    size = n_cells
  )
) +
  geom_point(alpha = 0.9, stroke = 0.4) +
  geom_text(
    aes(label = patient),
    size = 3,
    vjust = -0.7,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(3, 9)) +
  labs(
    title = "B. Patient switching landscape",
    x = "Within switching",
    y = "Between selection",
    color = "Histology",
    shape = "Condition",
    size = "Cell number"
  ) +
  theme_candle +
  theme(
    panel.grid.major = element_line(color = "grey92")
  )

# ----------------------------
# 7. Panel C:
#    Switch fraction by condition
# ----------------------------
pC <- ggplot(df, aes(x = condition, y = switch_fraction)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.2,
    alpha = 0.9
  ) +
  labs(
    title = "C. Switch fraction by condition",
    x = NULL,
    y = "Switch fraction"
  ) +
  theme_candle +
  theme(
    legend.position = "none"
  )

# ----------------------------
# 8. Panel D:
#    Switch fraction by histology
# ----------------------------
pD <- ggplot(df, aes(x = histology, y = switch_fraction, color = histology)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.2,
    alpha = 0.9
  ) +
  labs(
    title = "D. Switch fraction by histology",
    x = NULL,
    y = "Switch fraction"
  ) +
  theme_candle +
  theme(
    legend.position = "none"
  )

# ----------------------------
# 9. Combine into one panel
# ----------------------------
final_plot <-
  (pA | pB)+ #/
  #(pC | pD) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

# ----------------------------
# 10. Save outputs
# ----------------------------
ggsave(
  filename = "wilms_CANDLE_per_patient_summary_panel.pdf",
  plot = final_plot,
  width = 20,
  height = 7,
  units = "in",
  dpi = 300#,
 # device = cairo_pdf
)

# ggsave(
#   filename = "wilms_CANDLE_per_patient_summary_panel.png",
#   plot = final_plot,
#   width = 16,
#   height = 10,
#   units = "in",
#   dpi = 300
# )

# ----------------------------
# 11. Optional: print to screen
# ----------------------------
print(final_plot)

##### Sample-level CNV from numbat object

# obj <- load('integrated_compartment_cancer_processed.RData')
# obj <- s_objs
# colnames(obj@meta.data)
# table(obj@meta.data$subcompartment)
# blastema <- subset(obj,subset=subcompartment%in%c('blastema'))
# colnames(blastema@meta.data)

obj <- load('TME.RData')
obj <- tme_
colnames(obj@meta.data)
table(obj@meta.data$subcompartment)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# Read Numbat joint posterior table
jp <- read_tsv("./numbat_output/W944/joint_post_1.tsv", show_col_types = FALSE)

# Basic check
req <- c("cell", "CHROM", "seg", "cnv_state_map",
         "p_amp", "p_neu", "p_del", "p_loh", "p_bamp", "p_bdel",
         "p_cnv", "avg_entropy", "seg_start", "seg_end")
miss <- setdiff(req, colnames(jp))
if (length(miss) > 0) {
  stop("Missing columns: ", paste(miss, collapse = ", "))
}

# Optional: define segment length
jp <- jp %>%
  mutate(
    seg_len = seg_end - seg_start + 1,
    seg_cnv_prob = p_amp + p_del + p_loh + p_bamp + p_bdel,
    is_cnv = cnv_state_map != "neu",
    is_loh = cnv_state_map %in% c("loh", "bdel")
  )

# Optional signed state score
state_score_map <- c(
  "neu"  =  0,
  "amp"  =  1,
  "bamp" =  2,
  "del"  = -1,
  "bdel" = -2,
  "loh"  = -0.5
)

jp <- jp %>%
  mutate(
    cnv_state_score = unname(state_score_map[cnv_state_map])
  )

# Per-cell summary
cell_cnv <- jp %>%
  group_by(cell) %>%
  summarise(
    n_segments = n(),
    cnv_avg_prob = mean(seg_cnv_prob, na.rm = TRUE),
    cnv_avg_prob_w = weighted.mean(seg_cnv_prob, w = seg_len, na.rm = TRUE),
    
    frac_cnv_segments = mean(is_cnv, na.rm = TRUE),
    frac_cnv_segments_w = weighted.mean(as.numeric(is_cnv), w = seg_len, na.rm = TRUE),
    
    frac_loh_segments = mean(is_loh, na.rm = TRUE),
    frac_loh_segments_w = weighted.mean(as.numeric(is_loh), w = seg_len, na.rm = TRUE),
    
    cnv_avg_state = mean(abs(cnv_state_score), na.rm = TRUE),
    cnv_signed_state = mean(cnv_state_score, na.rm = TRUE),
    
    mean_p_cnv = mean(p_cnv, na.rm = TRUE),
    mean_entropy = mean(avg_entropy, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(cell_cnv, "./numbat_output/W944/joint_post_1.per_cell_cnv.tsv")

# suppressPackageStartupMessages({
#   library(Seurat)
#   library(dplyr)
#   library(readr)
#   library(stringr)
#   library(tibble)
#   library(ggplot2)
#   library(patchwork)
#   library(ComplexHeatmap)
#   library(circlize)
#   library(grid)
# })
# 
# # =========================================================
# # 1. LOAD OBJECT
# # =========================================================
# load("TME.RData")
# obj <- tme_
# 
# # =========================================================
# # 2. USER INPUTS
# # =========================================================
# sample_id <- "W001"
# numbat_root <- "./numbat_output"
# keep_subcompartments <- c("blastema", "epithelium", "stroma")
# state_col <- "subcompartment"
# reduction_use <- "umap"
# 
# preferred_cnv_cols <- c(
#   "cnv_avg_prob_w",
#   "cnv_avg_prob",
#   "frac_cnv_segments_w",
#   "frac_cnv_segments",
#   "cnv_avg_state",
#   "mean_p_cnv"
# )
# 
# # =========================================================
# # 3. HELPERS
# # =========================================================
# clean_barcode <- function(x) {
#   x %>%
#     as.character() %>%
#     str_replace("-1$", "") %>%
#     str_trim()
# }
# 
# extract_barcode_from_seurat_cell <- function(x) {
#   x <- as.character(x)
#   
#   x1 <- str_replace(x, "_[^_]+$", "")
#   x2 <- ifelse(str_detect(x1, "_"),
#                str_replace(x1, "^[^_]+_", ""),
#                x1)
#   
#   clean_barcode(x2)
# }
# 
# chr_to_factor <- function(x) {
#   x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
#   lev <- c(as.character(1:22), "X", "Y", "M", "MT")
#   factor(x, levels = lev, ordered = TRUE)
# }
# 
# # =========================================================
# # 4. SUBSET SAMPLE
# # =========================================================
# obj_s <- subset(obj, subset = orig.ident == sample_id)
# 
# cat("Cells in sample subset:", ncol(obj_s), "\n")
# print(table(obj_s$subcompartment, useNA = "ifany"))
# 
# # =========================================================
# # 5. READ PER-CELL CNV SUMMARY AND MERGE INTO SAMPLE SUBSET
# # =========================================================
# cnv_file <- file.path(numbat_root, sample_id, "joint_post_1.per_cell_cnv.tsv")
# if (!file.exists(cnv_file)) stop("Missing file: ", cnv_file)
# 
# cnv_df <- read_tsv(cnv_file, show_col_types = FALSE)
# if (!("cell" %in% colnames(cnv_df))) stop("Expected column 'cell' in ", cnv_file)
# 
# cnv_df <- cnv_df %>%
#   mutate(barcode_clean = clean_barcode(cell))
# 
# meta_s <- obj_s@meta.data %>%
#   rownames_to_column("seurat_cell") %>%
#   mutate(barcode_clean = extract_barcode_from_seurat_cell(seurat_cell))
# 
# cnv_cols_present <- intersect(preferred_cnv_cols, colnames(cnv_df))
# if (length(cnv_cols_present) == 0) {
#   stop("None of the preferred CNV summary columns found in ", cnv_file)
# }
# 
# meta_s2 <- meta_s %>%
#   left_join(
#     cnv_df %>% select(barcode_clean, all_of(cnv_cols_present)),
#     by = "barcode_clean"
#   )
# 
# rownames(meta_s2) <- meta_s2$seurat_cell
# obj_s@meta.data <- meta_s2[, setdiff(colnames(meta_s2), "seurat_cell"), drop = FALSE]
# 
# main_cnv_col <- cnv_cols_present[1]
# obj_s$cnv_burden_main <- obj_s@meta.data[[main_cnv_col]]
# 
# cat("Using main CNV column:", main_cnv_col, "\n")
# 
# # =========================================================
# # 6. RESTRICT TO BLASTEMA / EPITHELIUM / STROMA
# # =========================================================
# obj_fig <- subset(obj_s, subset = subcompartment %in% keep_subcompartments)
# 
# cat("Cells in figure subset:", ncol(obj_fig), "\n")
# print(table(obj_fig$subcompartment, useNA = "ifany"))
# 
# obj_fig_cnv <- subset(obj_fig, cells = colnames(obj_fig)[!is.na(obj_fig$cnv_burden_main)])
# 
# # =========================================================
# # 7. PANELS 1-3
# # =========================================================
# p1 <- DimPlot(
#   obj_fig,
#   reduction = reduction_use,
#   group.by = state_col,
#   label = TRUE,
#   repel = TRUE
# ) +
#   ggtitle(paste0(sample_id, ": cell states")) +
#   theme_bw(base_size = 12)
# 
# p2 <- FeaturePlot(
#   obj_fig_cnv,
#   reduction = reduction_use,
#   features = "cnv_burden_main"
# ) +
#   ggtitle(paste0(sample_id, ": CNV burden (", main_cnv_col, ")")) +
#   theme_bw(base_size = 12)
# 
# meta_plot <- obj_fig_cnv@meta.data %>%
#   rownames_to_column("cell") %>%
#   mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments))
# 
# p3 <- ggplot(meta_plot, aes(x = subcompartment, y = cnv_burden_main)) +
#   geom_violin(fill = "grey85", color = "black", trim = FALSE) +
#   geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
#   geom_jitter(width = 0.12, size = 0.35, alpha = 0.25) +
#   labs(
#     title = paste0(sample_id, ": CNV burden by subcompartment"),
#     x = NULL,
#     y = main_cnv_col
#   ) +
#   theme_bw(base_size = 12)
# 
# # =========================================================
# # 8. PANEL 4 = NUMBAT-LIKE CHROMOSOME CNV STATE HEATMAP
# # =========================================================
# joint_file <- file.path(numbat_root, sample_id, "joint_post_1.tsv")
# if (!file.exists(joint_file)) stop("Missing raw segment-level file: ", joint_file)
# 
# jp <- read_tsv(joint_file, show_col_types = FALSE)
# 
# req_cols <- c("cell", "CHROM", "seg", "seg_start", "seg_end", "cnv_state_map")
# miss <- setdiff(req_cols, colnames(jp))
# if (length(miss) > 0) {
#   stop("Missing required columns in joint_post_1.tsv: ", paste(miss, collapse = ", "))
# }
# 
# cells_keep <- obj_fig@meta.data %>%
#   rownames_to_column("seurat_cell") %>%
#   mutate(barcode_clean = extract_barcode_from_seurat_cell(seurat_cell)) %>%
#   select(seurat_cell, barcode_clean, subcompartment, cnv_burden_main)
# 
# jp2 <- jp %>%
#   mutate(barcode_clean = clean_barcode(cell)) %>%
#   inner_join(cells_keep, by = "barcode_clean")
# 
# seg_df <- jp2 %>%
#   distinct(CHROM, seg, seg_start, seg_end) %>%
#   mutate(
#     chr_fac = chr_to_factor(CHROM),
#     seg_id = paste0("chr", CHROM, ":", seg_start, "-", seg_end, "_seg", seg)
#   ) %>%
#   arrange(chr_fac, seg_start, seg_end)
# 
# jp2 <- jp2 %>%
#   mutate(
#     seg_id = paste0("chr", CHROM, ":", seg_start, "-", seg_end, "_seg", seg)
#   )
# 
# wide <- jp2 %>%
#   select(seg_id, seurat_cell, cnv_state_map) %>%
#   distinct() %>%
#   tidyr::pivot_wider(names_from = seurat_cell, values_from = cnv_state_map)
# 
# mat_df <- as.data.frame(wide)
# rownames(mat_df) <- mat_df$seg_id
# mat_df$seg_id <- NULL
# mat <- as.matrix(mat_df)
# 
# mat <- mat[seg_df$seg_id[seg_df$seg_id %in% rownames(mat)], , drop = FALSE]
# 
# cell_order_df <- cells_keep %>%
#   filter(seurat_cell %in% colnames(mat)) %>%
#   mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments)) %>%
#   arrange(subcompartment, desc(cnv_burden_main))
# 
# mat <- mat[, cell_order_df$seurat_cell, drop = FALSE]
# 
# state_colors <- c(
#   "blastema" = "#1b9e77",
#   "epithelium" = "#7570b3",
#   "stroma" = "#d95f02"
# )
# 
# ha <- HeatmapAnnotation(
#   subcompartment = cell_order_df$subcompartment,
#   col = list(subcompartment = state_colors),
#   show_annotation_name = TRUE
# )
# 
# row_chr <- seg_df %>%
#   filter(seg_id %in% rownames(mat)) %>%
#   arrange(match(seg_id, rownames(mat))) %>%
#   pull(CHROM)
# 
# row_chr <- factor(
#   gsub("^chr", "", row_chr, ignore.case = TRUE),
#   levels = c(as.character(1:22), "X", "Y", "M", "MT")
# )
# 
# cnv_colors <- c(
#   "bdel" = "#313695",
#   "del"  = "#74add1",
#   "loh"  = "#fdae61",
#   "neu"  = "#f7f7f7",
#   "amp"  = "#f46d43",
#   "bamp" = "#a50026"
# )
# 
# ht <- Heatmap(
#   mat,
#   name = "CNV state",
#   col = cnv_colors,
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   show_row_names = FALSE,
#   show_column_names = FALSE,
#   row_split = row_chr,
#   column_split = cell_order_df$subcompartment,
#   top_annotation = ha,
#   column_title = paste0(sample_id, ": segment-level CNV states"),
#   row_title = NULL,
#   use_raster = FALSE,
#   border = TRUE,
#   heatmap_legend_param = list(
#     at = names(cnv_colors),
#     labels = names(cnv_colors),
#     title = "CNV state"
#   )
# )
# 
# # =========================================================
# # 9. SAVE FIGURE WITHOUT grid.grabExpr()
# # =========================================================
# 
# out_prefix <- paste0("wilms_", sample_id, "_numbat_state_panel")
# 
# # Save top-left 3 ggplot panels as grobs
# g1 <- ggplotGrob(p1)
# g2 <- ggplotGrob(p2)
# g3 <- ggplotGrob(p3)
# 
# pdf(paste0(out_prefix, ".pdf"), width = 16, height = 12, onefile = TRUE)
# grid::grid.newpage()
# 
# # 2x2 layout
# lay <- grid::grid.layout(
#   nrow = 2, ncol = 2,
#   widths = unit(c(0.5, 0.5), "npc"),
#   heights = unit(c(0.45, 0.55), "npc")
# )
# 
# vp <- grid::viewport(layout = lay)
# grid::pushViewport(vp)
# 
# # panel 1
# grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
# grid::grid.draw(g1)
# grid::upViewport()
# 
# # panel 2
# grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
# grid::grid.draw(g2)
# grid::upViewport()
# 
# # panel 3
# grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
# grid::grid.draw(g3)
# grid::upViewport()
# 
# # panel 4 = heatmap drawn directly
# grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
# draw(ht, newpage = FALSE, merge_legends = TRUE)
# grid::upViewport(2)
# 
# dev.off()

# suppressPackageStartupMessages({
#   library(Seurat)
#   library(dplyr)
#   library(readr)
#   library(stringr)
#   library(tidyr)
#   library(tibble)
#   library(ComplexHeatmap)
#   library(circlize)
#   library(grid)
# })
# 
# # =========================================================
# # 1. LOAD OBJECT
# # =========================================================
# load("TME.RData")
# obj <- tme_
# 
# # =========================================================
# # 2. USER INPUTS
# # =========================================================
# sample_id <- "W001"   # change as needed
# numbat_root <- "./numbat_output"
# keep_subcompartments <- c("blastema", "epithelium", "stroma")
# 
# # output name
# out_pdf <- paste0("wilms_", sample_id, "_arm_level_CNV_heatmap.pdf")
# 
# # Optional: maximum cells per group for plotting
# # set to NULL to keep all
# max_cells_per_group <- 3000
# 
# # =========================================================
# # 3. HELPERS
# # =========================================================
# clean_barcode <- function(x) {
#   x %>%
#     as.character() %>%
#     str_replace("-1$", "") %>%
#     str_trim()
# }
# 
# extract_barcode_from_seurat_cell <- function(x) {
#   x <- as.character(x)
#   x1 <- str_replace(x, "_[^_]+$", "")
#   x2 <- ifelse(str_detect(x1, "_"),
#                str_replace(x1, "^[^_]+_", ""),
#                x1)
#   clean_barcode(x2)
# }
# 
# # Human centromere positions (approximate, hg38-compatible use)
# centromere_df <- tibble(
#   chrom = c(as.character(1:22), "X", "Y"),
#   centromere = c(
#     123400000, 93900000, 90900000, 50000000, 48800000, 59800000,
#     60100000, 45200000, 43000000, 39800000, 53400000, 35500000,
#     17700000, 17200000, 19000000, 36800000, 24000000, 17200000,
#     26500000, 27500000, 13200000, 14700000, 60600000, 12500000
#   )
# )
# 
# # order of chromosome arms
# arm_levels <- c(rbind(paste0(1:22, "p"), paste0(1:22, "q")))
# arm_levels <- c(arm_levels, "Xp", "Xq", "Yp", "Yq")
# 
# # state priority for hard collapsing
# state_priority <- c("bdel", "del", "loh", "neu", "amp", "bamp")
# 
# collapse_state_weighted <- function(state_vec, weight_vec) {
#   d <- tibble(state = state_vec, w = weight_vec) %>%
#     filter(!is.na(state), !is.na(w))
#   
#   if (nrow(d) == 0) return(NA_character_)
#   
#   s <- d %>%
#     group_by(state) %>%
#     summarise(w = sum(w), .groups = "drop")
#   
#   # choose state with max weighted coverage
#   s <- s %>%
#     mutate(priority = match(state, state_priority)) %>%
#     arrange(desc(w), priority)
#   
#   s$state[1]
# }
# 
# # =========================================================
# # 4. SUBSET OBJECT TO SAMPLE + TARGET SUBCOMPARTMENTS
# # =========================================================
# obj_s <- subset(
#   obj,
#   subset = orig.ident == sample_id & subcompartment %in% keep_subcompartments
# )
# 
# cat("Cells in subset:\n")
# print(table(obj_s$subcompartment, useNA = "ifany"))
# 
# meta_s <- obj_s@meta.data %>%
#   rownames_to_column("seurat_cell") %>%
#   mutate(
#     barcode_clean = extract_barcode_from_seurat_cell(seurat_cell),
#     subcompartment = factor(subcompartment, levels = keep_subcompartments)
#   ) %>%
#   select(seurat_cell, barcode_clean, subcompartment)
# 
# # optional downsampling for plotting
# if (!is.null(max_cells_per_group)) {
#   set.seed(1)
#   meta_s <- meta_s %>%
#     group_by(subcompartment) %>%
#     group_modify(~ {
#       n_keep <- min(nrow(.x), max_cells_per_group)
#       .x %>% slice_sample(n = n_keep)
#     }) %>%
#     ungroup()
# }
# 
# cat("Cells retained for heatmap:\n")
# print(table(meta_s$subcompartment))
# 
# # =========================================================
# # 5. READ RAW NUMBAT JOINT POSTERIOR
# # =========================================================
# joint_file <- file.path(numbat_root, sample_id, "joint_post_1.tsv")
# if (!file.exists(joint_file)) {
#   stop("Missing file: ", joint_file)
# }
# 
# jp <- read_tsv(joint_file, show_col_types = FALSE)
# 
# req_cols <- c("cell", "CHROM", "seg_start", "seg_end", "cnv_state_map")
# miss <- setdiff(req_cols, colnames(jp))
# if (length(miss) > 0) {
#   stop("Missing required columns: ", paste(miss, collapse = ", "))
# }
# 
# jp <- jp %>%
#   mutate(
#     barcode_clean = clean_barcode(cell),
#     chrom = gsub("^chr", "", as.character(CHROM), ignore.case = TRUE)
#   ) %>%
#   filter(chrom %in% c(as.character(1:22), "X", "Y")) %>%
#   inner_join(meta_s, by = "barcode_clean")
# 
# cat("Matched segment rows:", nrow(jp), "\n")
# 
# # =========================================================
# # 6. ASSIGN EACH SEGMENT TO p OR q ARM
# # =========================================================
# jp <- jp %>%
#   left_join(centromere_df, by = c("chrom")) %>%
#   mutate(
#     seg_mid = (seg_start + seg_end) / 2,
#     arm = ifelse(seg_mid < centromere, "p", "q"),
#     chrom_arm = paste0(chrom, arm),
#     seg_len = pmax(seg_end - seg_start + 1, 1)
#   )
# 
# jp$chrom_arm <- factor(jp$chrom_arm, levels = arm_levels, ordered = TRUE)
# 
# # =========================================================
# # 7. COLLAPSE SEGMENTS -> CHROMOSOME ARM STATE PER CELL
# # =========================================================
# arm_state_df <- jp %>%
#   filter(!is.na(chrom_arm)) %>%
#   mutate(chrom_arm = as.character(chrom_arm)) %>%
#   group_by(seurat_cell, subcompartment, chrom_arm) %>%
#   summarise(
#     arm_state = collapse_state_weighted(cnv_state_map, seg_len),
#     .groups = "drop"
#   )
# 
# # Build full cell x arm grid explicitly
# full_grid <- expand.grid(
#   seurat_cell = unique(meta_s$seurat_cell),
#   chrom_arm   = arm_levels,
#   stringsAsFactors = FALSE
# ) %>%
#   as_tibble() %>%
#   left_join(
#     meta_s %>% select(seurat_cell, subcompartment),
#     by = "seurat_cell"
#   ) %>%
#   left_join(
#     arm_state_df %>% mutate(chrom_arm = as.character(chrom_arm)),
#     by = c("seurat_cell", "subcompartment", "chrom_arm")
#   ) %>%
#   mutate(
#     arm_state = ifelse(is.na(arm_state), "neu", arm_state),
#     chrom_arm = factor(chrom_arm, levels = arm_levels, ordered = TRUE),
#     subcompartment = factor(subcompartment, levels = keep_subcompartments)
#   )
# 
# # =========================================================
# # 8. MAKE MATRIX: ROWS = CELLS, COLS = CHROMOSOME ARMS
# # =========================================================
# heat_df <- full_grid %>%
#   select(seurat_cell, subcompartment, chrom_arm, arm_state) %>%
#   tidyr::pivot_wider(
#     names_from = chrom_arm,
#     values_from = arm_state
#   )
# 
# row_annot <- heat_df %>%
#   select(seurat_cell, subcompartment)
# 
# mat_df <- heat_df %>%
#   select(-seurat_cell, -subcompartment)
# 
# mat <- as.matrix(mat_df)
# rownames(mat) <- row_annot$seurat_cell
# 
# # order rows by subcompartment
# row_order_df <- row_annot %>%
#   mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments)) %>%
#   arrange(subcompartment)
# 
# mat <- mat[row_order_df$seurat_cell, , drop = FALSE]
# 
# # =========================================================
# # 9. ANNOTATIONS AND COLORS
# # =========================================================
# subcompartment_colors <- c(
#   "blastema" = "#1b9e77",
#   "epithelium" = "#7570b3",
#   "stroma" = "#d95f02"
# )
# 
# cnv_colors <- c(
#   "bdel" = "#313695",
#   "del"  = "#74add1",
#   "loh"  = "#fdae61",
#   "neu"  = "#f7f7f7",
#   "amp"  = "#f46d43",
#   "bamp" = "#a50026"
# )
# 
# # top annotation with chromosome labels
# chr_labels <- gsub("(p|q)$", "", colnames(mat))
# arm_labels <- gsub("^[0-9XY]+", "", colnames(mat))
# 
# ha_top <- HeatmapAnnotation(
#   chromosome = anno_text(
#     chr_labels,
#     rot = 0,
#     gp = gpar(fontsize = 8),
#     just = "center",
#     location = 0.5
#   ),
#   arm = anno_text(
#     arm_labels,
#     rot = 0,
#     gp = gpar(fontsize = 7)
#   ),
#   annotation_name_side = "left",
#   show_annotation_name = FALSE
# )
# 
# ha_left <- rowAnnotation(
#   subcompartment = row_order_df$subcompartment,
#   col = list(subcompartment = subcompartment_colors),
#   show_annotation_name = FALSE
# )
# 
# # column splits by chromosome
# col_split <- factor(chr_labels, levels = c(as.character(1:22), "X", "Y"), ordered = TRUE)
# 
# # row split by subcompartment
# row_split <- factor(row_order_df$subcompartment, levels = keep_subcompartments)
# 
# # =========================================================
# # 10. DRAW HEATMAP
# # =========================================================
# ht <- Heatmap(
#   mat,
#   name = "CNV state",
#   col = cnv_colors,
#   rect_gp = gpar(col = NA),
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   show_row_names = FALSE,
#   show_column_names = FALSE,
#   row_split = row_split,
#   column_split = col_split,
#   left_annotation = ha_left,
#   top_annotation = ha_top,
#   row_title = "cells",
#   column_title = paste0(sample_id, ": chromosome-arm CNV states"),
#   border = TRUE,
#   use_raster = FALSE,
#   heatmap_legend_param = list(
#     at = names(cnv_colors),
#     labels = names(cnv_colors),
#     title = "CNV state"
#   )
# )
# 
# pdf(out_pdf, width = 18, height = 10)
# draw(
#   ht,
#   merge_legends = TRUE,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right"
# )
# dev.off()
# 
# cat("Saved:", out_pdf, "\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# =========================================================
# 1. LOAD OBJECT
# =========================================================
load("TME.RData")
obj <- tme_

# =========================================================
# 2. USER INPUTS
# =========================================================
sample_id <- "W001"
numbat_root <- "./numbat_output"

keep_subcompartments <- c("blastema", "epithelium", "stroma")
state_col <- "subcompartment"
reduction_use <- "umap"

# If NULL, keep all cells in each subcompartment block
# If e.g. 1000 or 2000, downsample each block for a cleaner heatmap
max_cells_per_group <- 1000

# preferred CNV summary columns from joint_post_1.per_cell_cnv.tsv
preferred_cnv_cols <- c(
  "cnv_avg_prob_w",
  "cnv_avg_prob",
  "frac_cnv_segments_w",
  "frac_cnv_segments",
  "cnv_avg_state",
  "mean_p_cnv"
)

# output prefix
out_prefix <- paste0("wilms_", sample_id, "_cnv")

# =========================================================
# 3. HELPERS
# =========================================================
clean_barcode <- function(x) {
  x %>%
    as.character() %>%
    str_replace("-1$", "") %>%
    str_trim()
}

# Try to recover the original barcode from Seurat cell names
# Handles common forms:
#   SAMPLE_AAAC...-1
#   AAAC...-1_SAMPLE
#   AAAC...-1
extract_barcode_from_seurat_cell <- function(x) {
  x <- as.character(x)
  
  x1 <- str_replace(x, "_[^_]+$", "")
  x2 <- ifelse(
    str_detect(x1, "_"),
    str_replace(x1, "^[^_]+_", ""),
    x1
  )
  
  clean_barcode(x2)
}

# approximate human centromere positions
centromere_df <- tibble(
  chrom = c(as.character(1:22), "X", "Y"),
  centromere = c(
    123400000, 93900000, 90900000, 50000000, 48800000, 59800000,
    60100000, 45200000, 43000000, 39800000, 53400000, 35500000,
    17700000, 17200000, 19000000, 36800000, 24000000, 17200000,
    26500000, 27500000, 13200000, 14700000, 60600000, 12500000
  )
)

# chromosome-arm ordering
arm_levels <- c(rbind(paste0(1:22, "p"), paste0(1:22, "q")))
arm_levels <- c(arm_levels, "Xp", "Xq", "Yp", "Yq")

# hard CNV state priority when weighted coverage ties are close
state_priority <- c("bdel", "del", "loh", "neu", "amp", "bamp")

collapse_state_weighted <- function(state_vec, weight_vec) {
  d <- tibble(state = state_vec, w = weight_vec) %>%
    filter(!is.na(state), !is.na(w))
  
  if (nrow(d) == 0) return(NA_character_)
  
  s <- d %>%
    group_by(state) %>%
    summarise(w = sum(w), .groups = "drop") %>%
    mutate(priority = match(state, state_priority)) %>%
    arrange(desc(w), priority)
  
  s$state[1]
}

# =========================================================
# 4. SUBSET TO ONE SAMPLE
# =========================================================
obj_s <- subset(obj, subset = orig.ident == sample_id)

cat("Cells in sample subset:", ncol(obj_s), "\n")
cat("Subcompartments in sample subset:\n")
print(table(obj_s$subcompartment, useNA = "ifany"))

# =========================================================
# 5. READ SAMPLE-SPECIFIC PER-CELL CNV SUMMARY
# =========================================================
cnv_file <- file.path(numbat_root, sample_id, "joint_post_1.per_cell_cnv.tsv")
if (!file.exists(cnv_file)) {
  stop("Missing file: ", cnv_file)
}

cnv_df <- read_tsv(cnv_file, show_col_types = FALSE)

if (!("cell" %in% colnames(cnv_df))) {
  stop("Expected column 'cell' in: ", cnv_file)
}

cnv_df <- cnv_df %>%
  mutate(barcode_clean = clean_barcode(cell))

cat("Columns in per-cell CNV file:\n")
print(colnames(cnv_df))

# =========================================================
# 6. MERGE PER-CELL CNV SUMMARY INTO SAMPLE SUBSET ONLY
# =========================================================
meta_s_all <- obj_s@meta.data %>%
  rownames_to_column("seurat_cell") %>%
  mutate(barcode_clean = extract_barcode_from_seurat_cell(seurat_cell))

cnv_cols_present <- intersect(preferred_cnv_cols, colnames(cnv_df))
if (length(cnv_cols_present) == 0) {
  stop("None of the preferred CNV summary columns were found in ", cnv_file)
}

meta_s_all2 <- meta_s_all %>%
  left_join(
    cnv_df %>% select(barcode_clean, all_of(cnv_cols_present)),
    by = "barcode_clean"
  )

rownames(meta_s_all2) <- meta_s_all2$seurat_cell
obj_s@meta.data <- meta_s_all2[, setdiff(colnames(meta_s_all2), "seurat_cell"), drop = FALSE]

main_cnv_col <- cnv_cols_present[1]
obj_s$cnv_burden_main <- obj_s@meta.data[[main_cnv_col]]

cat("\nUsing CNV burden column for figures: ", main_cnv_col, "\n", sep = "")
cat("Matched cells with non-missing CNV burden: ",
    sum(!is.na(obj_s$cnv_burden_main)), " / ", ncol(obj_s), "\n", sep = "")

# =========================================================
# 7. RESTRICT TO BLASTEMA / EPITHELIUM / STROMA
# =========================================================
obj_fig <- subset(obj_s, subset = subcompartment %in% keep_subcompartments)

cat("\nCells in figure subset:\n")
print(table(obj_fig$subcompartment, useNA = "ifany"))

obj_fig_cnv <- subset(
  obj_fig,
  cells = colnames(obj_fig)[!is.na(obj_fig$cnv_burden_main)]
)

cat("Cells with CNV burden available in figure subset:", ncol(obj_fig_cnv), "\n")

# metadata used for heatmap cell selection
meta_s <- obj_fig@meta.data %>%
  rownames_to_column("seurat_cell") %>%
  mutate(
    barcode_clean = extract_barcode_from_seurat_cell(seurat_cell),
    subcompartment = factor(subcompartment, levels = keep_subcompartments)
  ) %>%
  select(seurat_cell, barcode_clean, subcompartment)

# optional downsampling for heatmap readability
if (!is.null(max_cells_per_group)) {
  set.seed(1)
  meta_s <- meta_s %>%
    group_by(subcompartment) %>%
    group_modify(~ {
      n_keep <- min(nrow(.x), max_cells_per_group)
      dplyr::slice_sample(.x, n = n_keep)
    }) %>%
    ungroup()
}

cat("\nCells retained for heatmap after optional downsampling:\n")
print(table(meta_s$subcompartment))

# =========================================================
# 8. PANELS A-C
# =========================================================

# A. UMAP by subcompartment
p1 <- DimPlot(
  obj_fig,
  reduction = reduction_use,
  group.by = state_col,
  label = TRUE,
  repel = TRUE
) +
  ggtitle(paste0("A. ", sample_id, ": cell states")) +
  theme_bw(base_size = 12)

# B. UMAP by CNV burden
p2 <- FeaturePlot(
  obj_fig_cnv,
  reduction = reduction_use,
  features = "cnv_burden_main"
) +
  ggtitle(paste0("B. ", sample_id, ": CNV burden (", main_cnv_col, ")")) +
  theme_bw(base_size = 12)

# C. Violin plot by subcompartment
meta_plot <- obj_fig_cnv@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments))

p3 <- ggplot(meta_plot, aes(x = subcompartment, y = cnv_burden_main)) +
  geom_violin(fill = "grey85", color = "black", trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(width = 0.12, size = 0.35, alpha = 0.25) +
  labs(
    title = paste0("C. ", sample_id, ": CNV burden by subcompartment"),
    x = NULL,
    y = main_cnv_col
  ) +
  theme_bw(base_size = 12)

# =========================================================
# 9. PANEL D: CHROMOSOME-ARM HEATMAP FROM RAW joint_post_1.tsv
# =========================================================
joint_file <- file.path(numbat_root, sample_id, "joint_post_1.tsv")
if (!file.exists(joint_file)) {
  stop("Missing raw segment-level file: ", joint_file)
}

jp <- read_tsv(joint_file, show_col_types = FALSE)

req_cols <- c("cell", "CHROM", "seg_start", "seg_end", "cnv_state_map")
miss <- setdiff(req_cols, colnames(jp))
if (length(miss) > 0) {
  stop("Missing required columns in joint_post_1.tsv: ", paste(miss, collapse = ", "))
}

cells_keep <- meta_s %>%
  select(seurat_cell, barcode_clean, subcompartment)

jp <- jp %>%
  mutate(
    barcode_clean = clean_barcode(cell),
    chrom = gsub("^chr", "", as.character(CHROM), ignore.case = TRUE)
  ) %>%
  filter(chrom %in% c(as.character(1:22), "X", "Y")) %>%
  inner_join(cells_keep, by = "barcode_clean") %>%
  left_join(centromere_df, by = "chrom") %>%
  mutate(
    seg_mid = (seg_start + seg_end) / 2,
    arm = ifelse(seg_mid < centromere, "p", "q"),
    chrom_arm = paste0(chrom, arm),
    seg_len = pmax(seg_end - seg_start + 1, 1)
  )

cat("\nMatched raw segment rows for heatmap: ", nrow(jp), "\n", sep = "")

# collapse segments to chromosome-arm state per cell
arm_state_df <- jp %>%
  filter(!is.na(chrom_arm)) %>%
  mutate(chrom_arm = as.character(chrom_arm)) %>%
  group_by(seurat_cell, subcompartment, chrom_arm) %>%
  summarise(
    arm_state = collapse_state_weighted(cnv_state_map, seg_len),
    .groups = "drop"
  )

# explicit full cell x arm grid
full_grid <- expand.grid(
  seurat_cell = unique(cells_keep$seurat_cell),
  chrom_arm   = arm_levels,
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(
    cells_keep %>% select(seurat_cell, subcompartment),
    by = "seurat_cell"
  ) %>%
  left_join(
    arm_state_df %>% mutate(chrom_arm = as.character(chrom_arm)),
    by = c("seurat_cell", "subcompartment", "chrom_arm")
  ) %>%
  mutate(
    arm_state = ifelse(is.na(arm_state), "neu", arm_state),
    chrom_arm = factor(chrom_arm, levels = arm_levels, ordered = TRUE),
    subcompartment = factor(subcompartment, levels = keep_subcompartments)
  )

# make matrix: rows = cells, cols = chromosome arms
heat_df <- full_grid %>%
  select(seurat_cell, subcompartment, chrom_arm, arm_state) %>%
  pivot_wider(names_from = chrom_arm, values_from = arm_state)

row_annot <- heat_df %>%
  select(seurat_cell, subcompartment)

mat_df <- heat_df %>%
  select(-seurat_cell, -subcompartment)

mat <- as.matrix(mat_df)
rownames(mat) <- row_annot$seurat_cell

# row order by subcompartment
row_order_df <- row_annot %>%
  mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments)) %>%
  arrange(subcompartment)

mat <- mat[row_order_df$seurat_cell, , drop = FALSE]

# optional clustering within each block for cleaner structure
state_to_num <- c(
  "bdel" = -2,
  "del"  = -1,
  "loh"  = -0.5,
  "neu"  =  0,
  "amp"  =  1,
  "bamp" =  2
)

mat_num <- matrix(state_to_num[mat], nrow = nrow(mat), dimnames = dimnames(mat))

ord_within <- unlist(lapply(keep_subcompartments, function(sc) {
  idx <- which(row_order_df$subcompartment == sc)
  if (length(idx) <= 2) return(idx)
  idx[hclust(dist(mat_num[idx, , drop = FALSE]))$order]
}))

mat <- mat[ord_within, , drop = FALSE]
row_order_df <- row_order_df[ord_within, , drop = FALSE]

# annotations and colors
subcompartment_colors <- c(
  "blastema" = "#1b9e77",
  "epithelium" = "#7570b3",
  "stroma" = "#d95f02"
)

cnv_colors <- c(
  "bdel" = "#313695",
  "del"  = "#74add1",
  "loh"  = "#fdae61",
  "neu"  = "#f7f7f7",
  "amp"  = "#f46d43",
  "bamp" = "#a50026"
)

chr_labels <- gsub("(p|q)$", "", colnames(mat))
arm_labels <- gsub("^[0-9XY]+", "", colnames(mat))

ha_top <- HeatmapAnnotation(
  chr = anno_text(chr_labels, gp = gpar(fontsize = 7), just = "center"),
  arm = anno_text(arm_labels, gp = gpar(fontsize = 6), just = "center"),
  show_annotation_name = FALSE
)

ha_left <- rowAnnotation(
  subcompartment = row_order_df$subcompartment,
  col = list(subcompartment = subcompartment_colors),
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)

row_split <- factor(row_order_df$subcompartment, levels = keep_subcompartments)
col_split <- factor(chr_labels, levels = c(as.character(1:22), "X", "Y"), ordered = TRUE)

ht <- Heatmap(
  mat,
  name = "CNV state",
  col = cnv_colors,
  rect_gp = gpar(col = NA),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_split = row_split,
  column_split = col_split,
  left_annotation = ha_left,
  top_annotation = ha_top,
  row_title = NULL,
  column_title = paste0("D. ", sample_id, ": chromosome-arm CNV states"),
  border = TRUE,
  use_raster = FALSE,
  heatmap_legend_param = list(
    at = names(cnv_colors),
    labels = names(cnv_colors),
    title = "CNV state"
  )
)

# =========================================================
# PANEL E: AVERAGE CNV FROM EXACT PANEL D MATRIX
# IMPORTANT: enforce identical chromosome-arm column order
# =========================================================

mat_num_final <- matrix(
  state_to_num[mat],
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
)

avg_cnv_df <- as.data.frame(mat_num_final) %>%
  rownames_to_column("seurat_cell") %>%
  left_join(
    row_order_df %>% select(seurat_cell, subcompartment),
    by = "seurat_cell"
  ) %>%
  pivot_longer(
    cols = all_of(colnames(mat_num_final)),
    names_to = "chrom_arm",
    values_to = "cnv_numeric"
  ) %>%
  group_by(subcompartment, chrom_arm) %>%
  summarise(
    avg_cnv = mean(cnv_numeric, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    subcompartment = factor(subcompartment, levels = keep_subcompartments),
    chrom_arm = factor(chrom_arm, levels = colnames(mat_num_final), ordered = TRUE)
  )

avg_mat_df <- avg_cnv_df %>%
  select(subcompartment, chrom_arm, avg_cnv) %>%
  pivot_wider(
    names_from = chrom_arm,
    values_from = avg_cnv
  ) %>%
  arrange(subcompartment)

avg_row_names <- as.character(avg_mat_df$subcompartment)

avg_mat <- avg_mat_df %>%
  select(-subcompartment) %>%
  as.matrix()

rownames(avg_mat) <- avg_row_names

# CRITICAL FIX: force same column order as Panel D
avg_mat <- avg_mat[, colnames(mat_num_final), drop = FALSE]

# Make Panel E-specific labels from avg_mat itself
avg_chr_labels <- gsub("(p|q)$", "", colnames(avg_mat))
avg_arm_labels <- gsub("^[0-9XY]+", "", colnames(avg_mat))

ha_top_avg <- HeatmapAnnotation(
  chr = anno_text(avg_chr_labels, gp = gpar(fontsize = 7), just = "center"),
  arm = anno_text(avg_arm_labels, gp = gpar(fontsize = 6), just = "center"),
  show_annotation_name = FALSE
)

avg_col_split <- factor(
  avg_chr_labels,
  levels = c(as.character(1:22), "X", "Y"),
  ordered = TRUE
)

avg_cnv_col_fun <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#313695", "#74add1", "#f7f7f7", "#f46d43", "#a50026")
)

ht_avg <- Heatmap(
  avg_mat,
  name = "Mean CNV",
  col = avg_cnv_col_fun,
  rect_gp = gpar(col = "black", lwd = 0.25),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = avg_col_split,
  top_annotation = ha_top_avg,
  column_title = paste0("E. ", sample_id, ": average chromosome-arm CNV"),
  border = TRUE,
  heatmap_legend_param = list(
    title = "Mean CNV",
    at = c(-2, -1, 0, 1, 2),
    labels = c("bdel", "del", "neu", "amp", "bamp")
  )
)

# # =========================================================
# # 10. SAVE FINAL 4-PANEL FIGURE
# # =========================================================
# g1 <- ggplotGrob(p1)
# g2 <- ggplotGrob(p2)
# g3 <- ggplotGrob(p3)
# 
# out_pdf <- paste0(out_prefix, ".pdf")
# 
# pdf(out_pdf, width = 18, height = 12, onefile = TRUE)
# grid.newpage()
# 
# lay <- grid.layout(
#   nrow = 2, ncol = 2,
#   widths = unit(c(0.48, 0.52), "npc"),
#   heights = unit(c(0.42, 0.58), "npc")
# )
# 
# pushViewport(viewport(layout = lay))
# 
# # Panel A
# pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
# grid.draw(g1)
# upViewport()
# 
# # Panel B
# pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
# grid.draw(g2)
# upViewport()
# 
# # Panel C
# pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
# grid.draw(g3)
# upViewport()
# 
# # Panel D
# pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
# draw(
#   ht,
#   newpage = FALSE,
#   merge_legends = TRUE,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right"
# )
# upViewport(2)
# 
# dev.off()
# 
# cat("\nSaved 4-panel figure to: ", out_pdf, "\n", sep = "")

# =========================================================
# 10. SAVE FINAL 5-PANEL FIGURE
# =========================================================

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)

out_prefix <- paste0("wilms_", sample_id, "_cnv_panel")
out_pdf <- paste0(out_prefix, ".pdf")

pdf(out_pdf, width = 18, height = 15, onefile = TRUE)
grid.newpage()

lay <- grid.layout(
  nrow = 3, ncol = 2,
  widths = unit(c(0.48, 0.52), "npc"),
  heights = unit(c(0.32, 0.28, 0.40), "npc")
)

pushViewport(viewport(layout = lay))

# Panel A
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(g1)
upViewport()

# Panel B
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
grid.draw(g2)
upViewport()

# Panel C
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(g3)
upViewport()

# Panel E
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(
  ht_avg,
  newpage = FALSE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
upViewport()

# Panel D across bottom row
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
draw(
  ht,
  newpage = FALSE,
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
upViewport(2)

dev.off()

cat("\nSaved 5-panel figure to: ", out_pdf, "\n", sep = "")

# =========================================================
# 11. OPTIONAL: SAVE SUBSET OBJECTS / TABLES
# =========================================================
# saveRDS(obj_s, paste0(out_prefix, "_sample_subset_with_cnv.rds"))
# saveRDS(obj_fig, paste0(out_prefix, "_figure_subset.rds"))
# write_tsv(meta_plot, paste0(out_prefix, "_violin_metadata.tsv"))

cat("Saved sample subset object, figure subset object, and violin metadata table.\n")

###### With statistics
###### With statistics
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggpubr)
  library(patchwork)
})

# =========================================================
# 1. LOAD OBJECT
# =========================================================
load("TME.RData")
obj <- tme_

# =========================================================
# 2. USER INPUTS
# =========================================================
numbat_root <- "./numbat_output"

samples_use <- unique(obj$orig.ident)
# Or manually:
#samples_use <- c("W001")

keep_subcompartments <- c("blastema", "epithelium", "stroma")
state_col <- "subcompartment"
sample_col <- "orig.ident"
reduction_use <- "umap"

max_cells_per_group <- 1000

preferred_cnv_cols <- c(
  "cnv_avg_prob_w",
  "cnv_avg_prob",
  "frac_cnv_segments_w",
  "frac_cnv_segments",
  "cnv_avg_state",
  "mean_p_cnv"
)

outdir <- "wilms_CNV_panels_with_statistics"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 3. HELPERS
# =========================================================
clean_barcode <- function(x) {
  x %>%
    as.character() %>%
    str_replace("-1$", "") %>%
    str_trim()
}

extract_barcode_from_seurat_cell <- function(x) {
  x <- as.character(x)
  x1 <- str_replace(x, "_[^_]+$", "")
  x2 <- ifelse(
    str_detect(x1, "_"),
    str_replace(x1, "^[^_]+_", ""),
    x1
  )
  clean_barcode(x2)
}

centromere_df <- tibble(
  chrom = c(as.character(1:22), "X", "Y"),
  centromere = c(
    123400000, 93900000, 90900000, 50000000, 48800000, 59800000,
    60100000, 45200000, 43000000, 39800000, 53400000, 35500000,
    17700000, 17200000, 19000000, 36800000, 24000000, 17200000,
    26500000, 27500000, 13200000, 14700000, 60600000, 12500000
  )
)

arm_levels <- c(rbind(paste0(1:22, "p"), paste0(1:22, "q")))
arm_levels <- c(arm_levels, "Xp", "Xq", "Yp", "Yq")

state_priority <- c("bdel", "del", "loh", "neu", "amp", "bamp")

state_to_num <- c(
  "bdel" = -2,
  "del"  = -1,
  "loh"  = -0.5,
  "neu"  =  0,
  "amp"  =  1,
  "bamp" =  2
)

collapse_state_weighted <- function(state_vec, weight_vec) {
  d <- tibble(state = state_vec, w = weight_vec) %>%
    filter(!is.na(state), !is.na(w))
  
  if (nrow(d) == 0) return(NA_character_)
  
  s <- d %>%
    group_by(state) %>%
    summarise(w = sum(w), .groups = "drop") %>%
    mutate(priority = match(state, state_priority)) %>%
    arrange(desc(w), priority)
  
  s$state[1]
}

subcompartment_colors <- c(
  "blastema" = "#1b9e77",
  "epithelium" = "#7570b3",
  "stroma" = "#d95f02"
)

cnv_colors <- c(
  "bdel" = "#313695",
  "del"  = "#74add1",
  "loh"  = "#fdae61",
  "neu"  = "#f7f7f7",
  "amp"  = "#f46d43",
  "bamp" = "#a50026"
)

avg_cnv_col_fun <- circlize::colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("#313695", "#74add1", "#f7f7f7", "#f46d43", "#a50026")
)

# =========================================================
# 4. FUNCTION: PROCESS ONE SAMPLE
# =========================================================
process_one_sample <- function(sample_id) {
  
  message("\n==============================")
  message("Processing sample: ", sample_id)
  message("==============================")
  
  cnv_file <- file.path(numbat_root, sample_id, "joint_post_1.per_cell_cnv.tsv")
  joint_file <- file.path(numbat_root, sample_id, "joint_post_1.tsv")
  
  if (!file.exists(cnv_file)) {
    warning("Missing per-cell CNV file: ", cnv_file)
    return(NULL)
  }
  
  if (!file.exists(joint_file)) {
    warning("Missing joint CNV file: ", joint_file)
    return(NULL)
  }
  
  obj_s <- subset(obj, subset = orig.ident == sample_id)
  
  cnv_df <- read_tsv(cnv_file, show_col_types = FALSE)
  
  if (!("cell" %in% colnames(cnv_df))) {
    warning("No cell column in: ", cnv_file)
    return(NULL)
  }
  
  cnv_df <- cnv_df %>%
    mutate(barcode_clean = clean_barcode(cell))
  
  cnv_cols_present <- intersect(preferred_cnv_cols, colnames(cnv_df))
  
  if (length(cnv_cols_present) == 0) {
    warning("No preferred CNV summary column found for: ", sample_id)
    return(NULL)
  }
  
  main_cnv_col <- cnv_cols_present[1]
  
  meta_s_all <- obj_s@meta.data %>%
    rownames_to_column("seurat_cell") %>%
    mutate(barcode_clean = extract_barcode_from_seurat_cell(seurat_cell))
  
  meta_s_all2 <- meta_s_all %>%
    left_join(
      cnv_df %>% select(barcode_clean, all_of(cnv_cols_present)),
      by = "barcode_clean"
    )
  
  rownames(meta_s_all2) <- meta_s_all2$seurat_cell
  obj_s@meta.data <- meta_s_all2[, setdiff(colnames(meta_s_all2), "seurat_cell"), drop = FALSE]
  obj_s$cnv_burden_main <- obj_s@meta.data[[main_cnv_col]]
  
  obj_fig <- subset(obj_s, subset = subcompartment %in% keep_subcompartments)
  
  obj_fig_cnv <- subset(
    obj_fig,
    cells = colnames(obj_fig)[!is.na(obj_fig$cnv_burden_main)]
  )
  
  meta_plot <- obj_fig_cnv@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      sample = sample_id,
      subcompartment = factor(subcompartment, levels = keep_subcompartments)
    )
  
  # =======================================================
  # SAMPLE-WISE GLOBAL CNV BURDEN STATISTICS
  # =======================================================
  kw <- tryCatch({
    kruskal.test(cnv_burden_main ~ subcompartment, data = meta_plot)
  }, error = function(e) NULL)
  
  meta_stat <- meta_plot %>%
    filter(!is.na(cnv_burden_main), !is.na(subcompartment))
  
  pairwise_stats <- ggpubr::compare_means(
    formula = cnv_burden_main ~ subcompartment,
    data = meta_stat,
    method = "wilcox.test",
    p.adjust.method = "BH"
  ) %>%
    mutate(sample = sample_id)
  pairwise_stats <- pairwise_stats %>%
    mutate(
      p_label = paste0(p.signif, " (FDR=", signif(p.adj, 2), ")")
    )
  
  pairwise_stats_file <- file.path(outdir, paste0(sample_id, "_samplewise_pairwise_CNV_stats.tsv"))
  write_tsv(pairwise_stats, pairwise_stats_file)
  
  kw_label <- if (!is.null(kw)) {
    paste0("Kruskal-Wallis P = ", signif(kw$p.value, 3))
  } else {
    "Kruskal-Wallis P = NA"
  }
  
  # =======================================================
  # PANELS A-C
  # =======================================================
  p1 <- DimPlot(
    obj_fig,
    reduction = reduction_use,
    group.by = state_col,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0("A. ", sample_id, ": cell states")) +
    theme_bw(base_size = 12)
  
  p2 <- FeaturePlot(
    obj_fig_cnv,
    reduction = reduction_use,
    features = "cnv_burden_main"
  ) +
    ggtitle(paste0("B. ", sample_id, ": CNV burden (", main_cnv_col, ")")) +
    theme_bw(base_size = 12)
  
  p3 <- ggplot(meta_plot, aes(x = subcompartment, y = cnv_burden_main)) +
    geom_violin(fill = "grey85", color = "black", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(width = 0.12, size = 0.35, alpha = 0.25) +
    labs(
      title = paste0("C. ", sample_id, ": CNV burden by subcompartment"),
      subtitle = kw_label,
      x = NULL,
      y = main_cnv_col
    ) +
    theme_bw(base_size = 12)
  
  if (nrow(pairwise_stats) > 0 && length(unique(meta_stat$subcompartment)) >= 2) {
    p3 <- p3 +
      stat_compare_means(
        comparisons = list(
          c("blastema", "epithelium"),
          c("blastema", "stroma"),
          c("epithelium", "stroma")
        ),
        method = "wilcox.test",
        label = "p.format",
        hide.ns = FALSE
      )
  }
  
  # =======================================================
  # PANEL D MATRIX
  # =======================================================
  meta_s <- obj_fig@meta.data %>%
    rownames_to_column("seurat_cell") %>%
    mutate(
      barcode_clean = extract_barcode_from_seurat_cell(seurat_cell),
      subcompartment = factor(subcompartment, levels = keep_subcompartments)
    ) %>%
    select(seurat_cell, barcode_clean, subcompartment)
  
  if (!is.null(max_cells_per_group)) {
    set.seed(1)
    meta_s <- meta_s %>%
      group_by(subcompartment) %>%
      group_modify(~ {
        n_keep <- min(nrow(.x), max_cells_per_group)
        dplyr::slice_sample(.x, n = n_keep)
      }) %>%
      ungroup()
  }
  
  jp <- read_tsv(joint_file, show_col_types = FALSE)
  
  req_cols <- c("cell", "CHROM", "seg_start", "seg_end", "cnv_state_map")
  miss <- setdiff(req_cols, colnames(jp))
  
  if (length(miss) > 0) {
    warning("Missing required columns in joint_post_1.tsv for ", sample_id)
    return(NULL)
  }
  
  cells_keep <- meta_s %>%
    select(seurat_cell, barcode_clean, subcompartment)
  
  jp <- jp %>%
    mutate(
      barcode_clean = clean_barcode(cell),
      chrom = gsub("^chr", "", as.character(CHROM), ignore.case = TRUE)
    ) %>%
    filter(chrom %in% c(as.character(1:22), "X", "Y")) %>%
    inner_join(cells_keep, by = "barcode_clean") %>%
    left_join(centromere_df, by = "chrom") %>%
    mutate(
      seg_mid = (seg_start + seg_end) / 2,
      arm = ifelse(seg_mid < centromere, "p", "q"),
      chrom_arm = paste0(chrom, arm),
      seg_len = pmax(seg_end - seg_start + 1, 1)
    )
  
  arm_state_df <- jp %>%
    filter(!is.na(chrom_arm)) %>%
    mutate(chrom_arm = as.character(chrom_arm)) %>%
    group_by(seurat_cell, subcompartment, chrom_arm) %>%
    summarise(
      arm_state = collapse_state_weighted(cnv_state_map, seg_len),
      .groups = "drop"
    )
  
  full_grid <- expand.grid(
    seurat_cell = unique(cells_keep$seurat_cell),
    chrom_arm = arm_levels,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(
      cells_keep %>% select(seurat_cell, subcompartment),
      by = "seurat_cell"
    ) %>%
    left_join(
      arm_state_df %>% mutate(chrom_arm = as.character(chrom_arm)),
      by = c("seurat_cell", "subcompartment", "chrom_arm")
    ) %>%
    mutate(
      arm_state = ifelse(is.na(arm_state), "neu", arm_state),
      chrom_arm = factor(chrom_arm, levels = arm_levels, ordered = TRUE),
      subcompartment = factor(subcompartment, levels = keep_subcompartments)
    )
  
  heat_df <- full_grid %>%
    select(seurat_cell, subcompartment, chrom_arm, arm_state) %>%
    pivot_wider(names_from = chrom_arm, values_from = arm_state)
  
  row_annot <- heat_df %>%
    select(seurat_cell, subcompartment)
  
  mat_df <- heat_df %>%
    select(-seurat_cell, -subcompartment)
  
  mat <- as.matrix(mat_df)
  rownames(mat) <- row_annot$seurat_cell
  
  row_order_df <- row_annot %>%
    mutate(subcompartment = factor(subcompartment, levels = keep_subcompartments)) %>%
    arrange(subcompartment)
  
  mat <- mat[row_order_df$seurat_cell, , drop = FALSE]
  
  mat_num <- matrix(
    state_to_num[mat],
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  ord_within <- unlist(lapply(keep_subcompartments, function(sc) {
    idx <- which(row_order_df$subcompartment == sc)
    if (length(idx) <= 2) return(idx)
    idx[hclust(dist(mat_num[idx, , drop = FALSE]))$order]
  }))
  
  mat <- mat[ord_within, , drop = FALSE]
  row_order_df <- row_order_df[ord_within, , drop = FALSE]
  
  mat_num_final <- matrix(
    state_to_num[mat],
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  chr_labels <- gsub("(p|q)$", "", colnames(mat))
  arm_labels <- gsub("^[0-9XY]+", "", colnames(mat))
  
  ha_top <- HeatmapAnnotation(
    chr = anno_text(chr_labels, gp = gpar(fontsize = 7), just = "center"),
    arm = anno_text(arm_labels, gp = gpar(fontsize = 6), just = "center"),
    show_annotation_name = FALSE
  )
  
  ha_left <- rowAnnotation(
    subcompartment = row_order_df$subcompartment,
    col = list(subcompartment = subcompartment_colors),
    show_annotation_name = FALSE,
    width = unit(4, "mm")
  )
  
  row_split <- factor(row_order_df$subcompartment, levels = keep_subcompartments)
  col_split <- factor(chr_labels, levels = c(as.character(1:22), "X", "Y"), ordered = TRUE)
  
  ht <- Heatmap(
    mat,
    name = "CNV state",
    col = cnv_colors,
    rect_gp = gpar(col = NA),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = row_split,
    column_split = col_split,
    left_annotation = ha_left,
    top_annotation = ha_top,
    row_title = NULL,
    column_title = paste0("D. ", sample_id, ": chromosome-arm CNV states"),
    border = TRUE,
    use_raster = FALSE,
    heatmap_legend_param = list(
      at = names(cnv_colors),
      labels = names(cnv_colors),
      title = "CNV state"
    )
  )
  
  # =======================================================
  # PANEL E: AVERAGE CNV + SAMPLE-WISE STATISTICS
  # =======================================================
  avg_cnv_df <- as.data.frame(mat_num_final) %>%
    rownames_to_column("seurat_cell") %>%
    left_join(
      row_order_df %>% select(seurat_cell, subcompartment),
      by = "seurat_cell"
    ) %>%
    pivot_longer(
      cols = all_of(colnames(mat_num_final)),
      names_to = "chrom_arm",
      values_to = "cnv_numeric"
    ) %>%
    group_by(subcompartment, chrom_arm) %>%
    summarise(
      mean_cnv = mean(cnv_numeric, na.rm = TRUE),
      frac_gain = mean(cnv_numeric > 0, na.rm = TRUE),
      frac_loss = mean(cnv_numeric < 0, na.rm = TRUE),
      frac_altered = mean(cnv_numeric != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sample = sample_id,
      subcompartment = factor(subcompartment, levels = keep_subcompartments),
      chrom_arm = factor(chrom_arm, levels = colnames(mat_num_final), ordered = TRUE)
    )
  
  avg_cnv_file <- file.path(outdir, paste0(sample_id, "_panelE_average_CNV_by_state_arm.tsv"))
  write_tsv(avg_cnv_df, avg_cnv_file)
  
  avg_mat_df <- avg_cnv_df %>%
    select(subcompartment, chrom_arm, mean_cnv) %>%
    pivot_wider(names_from = chrom_arm, values_from = mean_cnv) %>%
    arrange(subcompartment)
  
  avg_row_names <- as.character(avg_mat_df$subcompartment)
  
  avg_mat <- avg_mat_df %>%
    select(-subcompartment) %>%
    as.matrix()
  
  rownames(avg_mat) <- avg_row_names
  avg_mat <- avg_mat[, colnames(mat_num_final), drop = FALSE]
  
  stopifnot(identical(colnames(avg_mat), colnames(mat_num_final)))
  
  # sample-wise arm-level Kruskal tests across states
  sample_arm_stats <- as.data.frame(mat_num_final) %>%
    rownames_to_column("seurat_cell") %>%
    left_join(
      row_order_df %>% select(seurat_cell, subcompartment),
      by = "seurat_cell"
    ) %>%
    pivot_longer(
      cols = all_of(colnames(mat_num_final)),
      names_to = "chrom_arm",
      values_to = "cnv_numeric"
    ) %>%
    group_by(chrom_arm) %>%
    summarise(
      p_kruskal = tryCatch(
        kruskal.test(cnv_numeric ~ subcompartment)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sample = sample_id,
      p_adj = p.adjust(p_kruskal, method = "BH"),
      signif = case_when(
        is.na(p_adj) ~ "",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      p_label = paste0(signif, " (FDR=", signif(p_adj, 2), ")")
    )
  
  write_tsv(
    sample_arm_stats,
    file.path(outdir, paste0(sample_id, "_samplewise_arm_level_CNV_stats.tsv"))
  )
  
  sig_labels <- sample_arm_stats %>%
    filter(signif != "") %>%
    mutate(
      x = match(chrom_arm, colnames(avg_mat)),
      label = signif
    )
  
  ht_avg <- Heatmap(
    avg_mat,
    name = "Mean CNV",
    col = avg_cnv_col_fun,
    rect_gp = gpar(col = "black", lwd = 0.25),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    column_split = col_split,
    top_annotation = ha_top,
    column_title = paste0("E. ", sample_id, ": average CNV by state; * FDR < 0.05"),
    border = TRUE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      arm <- colnames(avg_mat)[j]
      sig <- sample_arm_stats$signif[match(arm, sample_arm_stats$chrom_arm)]
      if (!is.na(sig) && sig != "" && i == 1) {
        grid.text(sig, x, y + unit(0.35, "npc"), gp = gpar(fontsize = 7, fontface = "bold"))
      }
    },
    heatmap_legend_param = list(
      title = "Mean CNV",
      at = c(-2, -1, 0, 1, 2),
      labels = c("bdel", "del", "neu", "amp", "bamp")
    )
  )
  
  # =======================================================
  # SAVE 5-PANEL PDF
  # =======================================================
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g3 <- ggplotGrob(p3)
  
  out_pdf <- file.path(outdir, paste0("wilms_", sample_id, "_collaborator_5panel_with_stats.pdf"))
  
  pdf(out_pdf, width = 18, height = 15, onefile = TRUE)
  grid.newpage()
  
  lay <- grid.layout(
    nrow = 3, ncol = 2,
    widths = unit(c(0.48, 0.52), "npc"),
    heights = unit(c(0.32, 0.28, 0.40), "npc")
  )
  
  pushViewport(viewport(layout = lay))
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(g1)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.draw(g2)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.draw(g3)
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  draw(
    ht_avg,
    newpage = FALSE,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  upViewport()
  
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:2))
  draw(
    ht,
    newpage = FALSE,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  upViewport(2)
  
  dev.off()
  
  message("Saved: ", out_pdf)
  
  return(list(
    sample = sample_id,
    meta_plot = meta_plot,
    avg_cnv_df = avg_cnv_df,
    sample_arm_stats = sample_arm_stats,
    pairwise_stats = pairwise_stats
  ))
}

# =========================================================
# 5. RUN ALL SAMPLES
# =========================================================
res_list <- lapply(samples_use, process_one_sample)
res_list <- Filter(Negate(is.null), res_list)

all_meta <- bind_rows(lapply(res_list, `[[`, "meta_plot"))
all_avg_cnv <- bind_rows(lapply(res_list, `[[`, "avg_cnv_df"))
all_sample_arm_stats <- bind_rows(lapply(res_list, `[[`, "sample_arm_stats"))
all_pairwise_stats <- bind_rows(lapply(res_list, `[[`, "pairwise_stats"))

write_tsv(all_meta, file.path(outdir, "ALL_samples_cell_level_CNV_metadata.tsv"))
write_tsv(all_avg_cnv, file.path(outdir, "ALL_samples_average_CNV_by_state_arm.tsv"))
write_tsv(all_sample_arm_stats, file.path(outdir, "ALL_samples_samplewise_arm_level_stats.tsv"))
write_tsv(all_pairwise_stats, file.path(outdir, "ALL_samples_pairwise_CNV_burden_stats.tsv"))

# # =========================================================
# # 6. GLOBAL STATISTICS
# # =========================================================
# 
# # Global per-sample mean CNV burden
# sample_state_summary <- all_meta %>%
#   filter(!is.na(cnv_burden_main), subcompartment %in% keep_subcompartments) %>%
#   group_by(sample, subcompartment) %>%
#   summarise(
#     mean_cnv_burden = mean(cnv_burden_main, na.rm = TRUE),
#     median_cnv_burden = median(cnv_burden_main, na.rm = TRUE),
#     n_cells = n(),
#     .groups = "drop"
#   )
# 
# write_tsv(sample_state_summary, file.path(outdir, "GLOBAL_sample_state_CNV_burden_summary.tsv"))
# 
# global_kw <- kruskal.test(mean_cnv_burden ~ subcompartment, data = sample_state_summary)
# 
# global_pairwise <- ggpubr::compare_means(
#   mean_cnv_burden ~ subcompartment,
#   data = sample_state_summary,
#   method = "wilcox.test",
#   p.adjust.method = "BH"
# )
# 
# write_tsv(global_pairwise, file.path(outdir, "GLOBAL_pairwise_CNV_burden_stats.tsv"))
# 
# # Global arm-level average CNV statistics using sample-level values
# global_arm_stats <- all_avg_cnv %>%
#   group_by(chrom_arm) %>%
#   summarise(
#     p_kruskal = tryCatch(
#       kruskal.test(mean_cnv ~ subcompartment)$p.value,
#       error = function(e) NA_real_
#     ),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     p_adj = p.adjust(p_kruskal, method = "BH"),
#     signif = case_when(
#       is.na(p_adj) ~ "",
#       p_adj < 0.001 ~ "***",
#       p_adj < 0.01 ~ "**",
#       p_adj < 0.05 ~ "*",
#       TRUE ~ ""
#     )
#   )
# 
# write_tsv(global_arm_stats, file.path(outdir, "GLOBAL_arm_level_CNV_stats.tsv"))
# 
# # Global average CNV matrix
# global_avg_mat_df <- all_avg_cnv %>%
#   group_by(subcompartment, chrom_arm) %>%
#   summarise(
#     mean_cnv = mean(mean_cnv, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     subcompartment = factor(subcompartment, levels = keep_subcompartments),
#     chrom_arm = factor(chrom_arm, levels = arm_levels, ordered = TRUE)
#   ) %>%
#   select(subcompartment, chrom_arm, mean_cnv) %>%
#   pivot_wider(names_from = chrom_arm, values_from = mean_cnv) %>%
#   arrange(subcompartment)
# 
# global_row_names <- as.character(global_avg_mat_df$subcompartment)
# 
# global_avg_mat <- global_avg_mat_df %>%
#   select(-subcompartment) %>%
#   as.matrix()
# 
# rownames(global_avg_mat) <- global_row_names
# global_avg_mat <- global_avg_mat[, arm_levels, drop = FALSE]
# 
# global_chr_labels <- gsub("(p|q)$", "", colnames(global_avg_mat))
# global_arm_labels <- gsub("^[0-9XY]+", "", colnames(global_avg_mat))
# 
# global_ha_top <- HeatmapAnnotation(
#   chr = anno_text(global_chr_labels, gp = gpar(fontsize = 7), just = "center"),
#   arm = anno_text(global_arm_labels, gp = gpar(fontsize = 6), just = "center"),
#   show_annotation_name = FALSE
# )
# 
# global_col_split <- factor(
#   global_chr_labels,
#   levels = c(as.character(1:22), "X", "Y"),
#   ordered = TRUE
# )
# 
# ht_global <- Heatmap(
#   global_avg_mat,
#   name = "Global mean CNV",
#   col = avg_cnv_col_fun,
#   rect_gp = gpar(col = "black", lwd = 0.25),
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   show_row_names = TRUE,
#   show_column_names = FALSE,
#   column_split = global_col_split,
#   top_annotation = global_ha_top,
#   column_title = "Global average chromosome-arm CNV by cell state; * FDR < 0.05",
#   border = TRUE,
#   cell_fun = function(j, i, x, y, width, height, fill) {
#     arm <- colnames(global_avg_mat)[j]
#     sig <- global_arm_stats$signif[match(arm, global_arm_stats$chrom_arm)]
#     if (!is.na(sig) && sig != "" && i == 1) {
#       grid.text(sig, x, y + unit(0.35, "npc"), gp = gpar(fontsize = 7, fontface = "bold"))
#     }
#   },
#   heatmap_legend_param = list(
#     title = "Global mean CNV",
#     at = c(-2, -1, 0, 1, 2),
#     labels = c("bdel", "del", "neu", "amp", "bamp")
#   )
# )
# 
# # =========================================================
# # 7. GLOBAL PDF
# # =========================================================
# 
# p_global_burden <- ggplot(
#   sample_state_summary,
#   aes(x = subcompartment, y = mean_cnv_burden)
# ) +
#   geom_boxplot(outlier.shape = NA, fill = "grey85", color = "black") +
#   geom_jitter(aes(shape = sample), width = 0.12, size = 2.5, alpha = 0.8) +
#   stat_compare_means(
#     method = "kruskal.test",
#     label.y = max(sample_state_summary$mean_cnv_burden, na.rm = TRUE) * 1.1
#   ) +
#   stat_compare_means(
#     comparisons = list(
#       c("blastema", "epithelium"),
#       c("blastema", "stroma"),
#       c("epithelium", "stroma")
#     ),
#     method = "wilcox.test",
#     label = "p.signif",
#     hide.ns = FALSE
#   ) +
#   labs(
#     title = "Global CNV burden by cell state",
#     subtitle = paste0("Sample-level means; Kruskal-Wallis P = ", signif(global_kw$p.value, 3)),
#     x = NULL,
#     y = "Mean CNV burden per sample"
#   ) +
#   theme_bw(base_size = 12)
# 
# global_pdf <- file.path(outdir, "GLOBAL_CNV_statistics.pdf")
# 
# pdf(global_pdf, width = 18, height = 10, onefile = TRUE)
# grid.newpage()
# 
# lay <- grid.layout(
#   nrow = 1,
#   ncol = 2,
#   widths = unit(c(0.38, 0.62), "npc")
# )
# 
# pushViewport(viewport(layout = lay))
# 
# pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
# grid.draw(ggplotGrob(p_global_burden))
# upViewport()
# 
# pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
# draw(
#   ht_global,
#   newpage = FALSE,
#   merge_legends = TRUE,
#   heatmap_legend_side = "right",
#   annotation_legend_side = "right"
# )
# upViewport(2)
# 
# dev.off()
# 
# message("\nSaved global statistics PDF: ", global_pdf)
# message("All outputs saved in: ", outdir)
# 
