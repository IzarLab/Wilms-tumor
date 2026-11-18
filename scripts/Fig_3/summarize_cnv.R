#!/usr/bin/env Rscript

# Summarize first-iteration Numbat posteriors.
# Please change the path.
#
# Typical Usage from the repository root
#   Rscript scripts/Fig_3/summarize_cnv.R \
#     data/compartments/TME.RData numbat_output out/cnv
#
# Inputs, in positional order
#   1. TME.RData: a Seurat object named tme_ with nucleus names, orig.ident
#      (sample ID), subcompartment (cell annotation) and a umap reduction
#      for the per-sample diagnostic plots.
#   2. NUMBAT_ROOT: one subdirectory per sample, named by orig.ident.
#      Each sample requires joint_post_1.tsv (first-iteration posteriors).
#      If supplied, joint_post_1.per_cell_cnv.tsv is reused. Otherwise the
#      per-nucleus burden table is computed from the joint posteriors.
#   3. OUTPUT_DIR: destination for cohort-level tables and diagnostic plots.
#
# Required R packages
#   Seurat, dplyr, readr, stringr, tidyr, tibble, ggplot2, ComplexHeatmap,
#   circlize and ggpubr. Install ComplexHeatmap via Bioconductor.
#
# Outputs consumed by scripts/Fig_3/a.R
#   ALL_samples_average_CNV_by_state_arm.tsv
#   ALL_samples_samplewise_arm_level_stats.tsv
#   ALL_samples_cell_level_CNV_metadata.tsv
#   The figure script reads these three tables from data/Misc/.
#   This script also writes pairwise burden statistics and sample PDFs.
#
# Analysis scope
#   The posterior files can contain malignant and nonmalignant nuclei.
#   The figure summaries use annotated blastemal, epithelial and stromal
#   nuclei. Burden plots use nuclei with available burden values. Arm
#   summaries use up to 1,000 nuclei per compartment per sample (seed 1).
#   Missing arm calls use the original neutral-fill convention.
#   The segment coordinates must be compatible with the centromere
#   coordinates below before applying the arm classification.
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
})

# =========================================================
# 1. LOAD OBJECT
# =========================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/Fig_3/cnv/summarize_cnv.R TME.RData NUMBAT_ROOT OUTPUT_DIR")
}
input_data <- new.env()
load(args[1], envir = input_data)
obj <- input_data$tme_
if (!inherits(obj, "Seurat") ||
    !all(c("orig.ident", "subcompartment") %in% names(obj@meta.data)) ||
    !"umap" %in% names(obj@reductions)) {
  stop("TME.RData must contain tme_ with orig.ident, subcompartment and a umap reduction.")
}

# =========================================================
# 2. USER INPUTS
# =========================================================
numbat_root <- args[2]

samples_use <- unique(obj$orig.ident)
# To run selected samples, replace samples_use with their orig.ident values.

keep_subcompartments <- c("blastema", "epithelium", "stroma")
state_col <- "subcompartment"
reduction_use <- "umap"

max_cells_per_group <- 1000

# Historical burden tables may contain only alternative summary columns.
# Use cnv_avg_prob_w whenever it is available and report any fallback.
preferred_cnv_cols <- c(
  "cnv_avg_prob_w",
  "cnv_avg_prob",
  "frac_cnv_segments_w",
  "frac_cnv_segments",
  "cnv_avg_state",
  "mean_p_cnv"
)

outdir <- args[3]
if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE)) {
  stop("Could not create output directory: ", outdir)
}

# =========================================================
# 3. HELPERS
# =========================================================
clean_barcode <- function(x) {
  x %>%
    as.character() %>%
    str_replace("-1$", "") %>%
    str_trim()
}

# Remove the Wilms Seurat naming wrappers before comparing with the
# per-sample Numbat barcodes. The join is scoped to one sample at a time.
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

# Approximate hg38 centromere positions from the original Wilms script.
# Include X and Y on the arm grid alongside the autosomes.
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

# Break ties in total segment length using this state order.
state_priority <- c("bdel", "del", "loh", "neu", "amp", "bamp")

# Signed plotting scores, expressed on the original ordinal scale.
# Copy-neutral LOH carries a score of -0.5 to distinguish it from neutrality.
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
summarize_burden <- function(jp) {
  # Check the posterior and coordinate fields used by the burden summaries.
  req <- c("cell", "CHROM", "seg", "cnv_state_map",
           "p_amp", "p_neu", "p_del", "p_loh", "p_bamp", "p_bdel",
           "p_cnv", "avg_entropy", "seg_start", "seg_end")
  miss <- setdiff(req, colnames(jp))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }

  # Sum the probabilities of all five non-neutral states, including
  # copy-neutral LOH. Use the inclusive segment length in base pairs.
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

  # cnv_avg_prob_w is the weighted mean of the non-neutral posterior sum.
  # Segments with missing posterior sums are omitted from that mean.
  # The other columns preserve the historical diagnostic summaries.
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
  cell_cnv
}

process_one_sample <- function(sample_id) {

  message("\n==============================")
  message("Processing sample: ", sample_id)
  message("==============================")

  cnv_file <- file.path(numbat_root, sample_id, "joint_post_1.per_cell_cnv.tsv")
  joint_file <- file.path(numbat_root, sample_id, "joint_post_1.tsv")

  if (!file.exists(joint_file)) {
    warning("Missing joint CNV file: ", joint_file)
    return(NULL)
  }

  obj_s <- subset(obj, subset = orig.ident == sample_id)

  if (file.exists(cnv_file)) {
    cnv_df <- read_tsv(cnv_file, show_col_types = FALSE)
  } else {
    cnv_df <- summarize_burden(read_tsv(joint_file, show_col_types = FALSE))
    write_tsv(cnv_df, file.path(outdir, paste0(sample_id, "_per_cell_cnv.tsv")))
  }

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
  if (main_cnv_col != "cnv_avg_prob_w") {
    warning("Using ", main_cnv_col, " as the burden column for ", sample_id)
  }

  meta_s_all <- obj_s@meta.data %>%
    rownames_to_column("seurat_cell") %>%
    mutate(barcode_clean = extract_barcode_from_seurat_cell(seurat_cell))

  if (anyDuplicated(meta_s_all$barcode_clean) || anyDuplicated(cnv_df$barcode_clean)) {
    stop("Duplicate barcodes after sample-specific normalization for ", sample_id)
  }

  # Preserve annotation rows even when posterior-based burden is unavailable.
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

  # Downsample the arm comparison only. Use the same seed for each sample
  # and preserve the original order of compartment processing.
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

  # Assign each segment entirely to the arm containing its midpoint.
  # A midpoint at the centromere is assigned to the q arm.
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

  if (!nrow(jp)) {
    stop("No posterior segment rows matched sampled malignant barcodes for ", sample_id)
  }

  # Pick one MAP state per nucleus and arm by summing lengths by state.
  arm_state_df <- jp %>%
    filter(!is.na(chrom_arm)) %>%
    mutate(chrom_arm = as.character(chrom_arm)) %>%
    group_by(seurat_cell, subcompartment, chrom_arm) %>%
    summarise(
      arm_state = collapse_state_weighted(cnv_state_map, seg_len),
      .groups = "drop"
    )

  # Expand to every sampled nucleus and arm. Arms without a matched state
  # receive a neutral score, including sex-chromosome arms without calls.
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
  # Each sampled nucleus contributes one score per arm. Segment length
  # determines its categorical arm call. The group mean weights nuclei
  # equally, and the signed fraction summaries count LOH as negative.
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

  # Compare component-level arm score distributions within each tumor.
  # Adjust the Kruskal-Wallis tests across arms in that tumor by BH.
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
# Sample directories without joint posteriors are reported and skipped.
# Cohort tables contain only the samples that could be summarized.
res_list <- lapply(samples_use, process_one_sample)
res_list <- Filter(Negate(is.null), res_list)
if (!length(res_list)) stop("No samples were summarized.")

all_meta <- bind_rows(lapply(res_list, `[[`, "meta_plot"))
all_avg_cnv <- bind_rows(lapply(res_list, `[[`, "avg_cnv_df"))
all_sample_arm_stats <- bind_rows(lapply(res_list, `[[`, "sample_arm_stats"))
all_pairwise_stats <- bind_rows(lapply(res_list, `[[`, "pairwise_stats"))

write_tsv(all_meta, file.path(outdir, "ALL_samples_cell_level_CNV_metadata.tsv"))
write_tsv(all_avg_cnv, file.path(outdir, "ALL_samples_average_CNV_by_state_arm.tsv"))
write_tsv(all_sample_arm_stats, file.path(outdir, "ALL_samples_samplewise_arm_level_stats.tsv"))
write_tsv(all_pairwise_stats, file.path(outdir, "ALL_samples_pairwise_CNV_burden_stats.tsv"))
