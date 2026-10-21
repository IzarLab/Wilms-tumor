# ============================================================
# Treatment-regimen-specific malignant trajectory analysis
# ============================================================
#
# GitHub script name : 02_analyze_treatment_regimen_specific_monocle_trajectories.R
# Original file      : Fig.3.monocle.wilms.treatment.regimes.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Runs Monocle3 trajectory/transition analysis across EE-4A, DD-4A, and Other treated regimens.
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
#   results/figure_04/
#
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(pheatmap)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

set.seed(1234)

# ============================================================
# 0. SETTINGS
# ============================================================

obj_path <- "/home/ubuntu/wilms/integrated_compartment_cancer_processed.RData"

outdir <- "./Wilms_treated_regimen_EE4A_DD4A_Other_shared_trajectory"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

condition_col <- "condition"
sample_col    <- "orig.ident"
celltype_col  <- "cell_type"

treated_condition <- "neoadjuvant"

assay_counts   <- "RNA"
umap_reduction <- "umap"

root_celltype <- "blastema 1"

regimen_map <- c(
  "W002" = "Other",
  "W010" = "EE-4A",
  "W009" = "EE-4A",
  "W003" = "EE-4A",
  "W441" = "EE-4A",
  "W944" = "DD-4A",
  "W184" = "DD-4A",
  "W180" = "DD-4A",
  "W370" = "DD-4A",
  "W006" = "DD-4A",
  "W050" = "Other",
  "W001" = "Other"
)

regimen_levels <- c("EE-4A", "DD-4A", "Other")

regimen_colors <- c(
  "EE-4A" = "#D73027",
  "DD-4A" = "#1A9850",
  "Other" = "#8C510A"
)

celltype_order <- c(
  "blastema 1",
  "S blastema 1",
  "G2M blastema 1",
  "blastema 2",
  "G2M blastema 2",
  "blastema 3",
  "smooth myocyte-like precursors",
  "G2M smooth myocyte-like precursors",
  "nascent smooth myocyte-like",
  "smooth myocyte-like",
  "S smooth myocyte-like",
  "G2M smooth myocyte-like",
  "PAX3+ myogenic precursors",
  "myoblast-like",
  "striated myocyte-like",
  "tubules",
  "S tubules",
  "G2M tubules",
  "podocyte-like"
)

n_per_celltype_per_regimen <- 1200

gene_min_cells     <- 10
cluster_k          <- 20
cluster_resolution <- 1e-3
minimal_branch_len <- 10

k_neighbors <- 20
pseudotime_forward_min_delta <- 0.01
max_neighbors_used <- 8

min_flux_to_draw <- 0.005

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

msg <- function(...) cat(sprintf(...), "\n")

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(NA_real_, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

get_counts_safe <- function(obj, assay = "RNA") {
  assays_available <- names(obj@assays)
  
  if (!assay %in% assays_available) {
    stop(sprintf(
      "Assay '%s' not found. Available assays: %s",
      assay,
      paste(assays_available, collapse = ", ")
    ))
  }
  
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  
  if (is.null(mat)) {
    mat <- tryCatch(
      GetAssayData(obj, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  
  if (is.null(mat)) {
    stop(sprintf("Could not extract counts from assay '%s'.", assay))
  }
  
  mat
}

save_gg <- function(filename, plot, width = 8, height = 6) {
  ggsave(
    filename,
    plot = plot,
    width = width,
    height = height,
    dpi = 300,
    device = cairo_pdf
  )
}

get_knn_indices <- function(mat, k = 20) {
  d <- as.matrix(dist(mat))
  diag(d) <- Inf
  idx <- t(apply(d, 1, function(z) order(z)[seq_len(min(k, length(z)))]))
  idx
}

matrix_to_long <- function(mat, value_name = "value") {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    rename(source = Var1, target = Var2, !!value_name := Freq) %>%
    mutate(
      source = as.character(source),
      target = as.character(target)
    )
}

build_cds_with_seurat_umap <- function(
    obj_small,
    assay_counts = "RNA",
    umap_reduction = "umap",
    gene_min_cells = 10,
    cluster_k = 20,
    cluster_resolution = 1e-3,
    minimal_branch_len = 10
) {
  counts_mat <- get_counts_safe(obj_small, assay = assay_counts)
  
  cell_metadata <- obj_small@meta.data
  
  gene_metadata <- data.frame(
    gene_short_name = rownames(counts_mat),
    row.names = rownames(counts_mat)
  )
  
  cds <- new_cell_data_set(
    expression_data = counts_mat,
    cell_metadata   = cell_metadata,
    gene_metadata   = gene_metadata
  )
  
  cds <- cds[rowSums(counts(cds) > 0) >= gene_min_cells, ]
  
  umap_mat <- Embeddings(obj_small, umap_reduction)
  umap_mat <- umap_mat[colnames(obj_small), , drop = FALSE]
  colnames(umap_mat)[1:2] <- c("UMAP_1", "UMAP_2")
  
  reducedDims(cds)$UMAP <- umap_mat
  
  cds <- cluster_cells(
    cds,
    reduction_method = "UMAP",
    k = cluster_k,
    resolution = cluster_resolution
  )
  
  cds <- learn_graph(
    cds,
    use_partition = FALSE,
    close_loop = FALSE,
    learn_graph_control = list(minimal_branch_len = minimal_branch_len)
  )
  
  cds
}

compute_directed_transition_flux <- function(
    df_condition,
    subtype_col = "cell_type",
    pt_col = "pseudotime",
    k_neighbors = 20,
    min_delta = 0.01,
    max_neighbors_used = 8,
    state_levels = NULL
) {
  df_condition <- df_condition %>%
    filter(!is.na(.data[[subtype_col]]), !is.na(.data[[pt_col]])) %>%
    distinct(cell, .keep_all = TRUE)
  
  if (nrow(df_condition) < 5) {
    return(list(edge_table = data.frame(), flux_matrix = matrix(0, 0, 0)))
  }
  
  if (is.null(state_levels)) {
    state_levels <- sort(unique(df_condition[[subtype_col]]))
  }
  
  coords <- as.matrix(df_condition[, c("UMAP_1", "UMAP_2")])
  knn_idx <- get_knn_indices(coords, k = min(k_neighbors, nrow(df_condition) - 1))
  
  edge_list <- vector("list", nrow(df_condition))
  
  for (i in seq_len(nrow(df_condition))) {
    nbr_idx <- knn_idx[i, ]
    nbr_idx <- nbr_idx[!is.na(nbr_idx)]
    
    pt_i <- df_condition[[pt_col]][i]
    pt_j <- df_condition[[pt_col]][nbr_idx]
    
    keep <- which((pt_j - pt_i) >= min_delta)
    
    if (length(keep) == 0) {
      edge_list[[i]] <- NULL
      next
    }
    
    nbr_idx <- nbr_idx[keep]
    pt_j <- pt_j[keep]
    
    ord <- order(pt_j - pt_i)
    nbr_idx <- nbr_idx[ord][seq_len(min(length(ord), max_neighbors_used))]
    
    w <- rep(1 / length(nbr_idx), length(nbr_idx))
    
    edge_list[[i]] <- data.frame(
      source_cell  = df_condition$cell[i],
      target_cell  = df_condition$cell[nbr_idx],
      source_state = df_condition[[subtype_col]][i],
      target_state = df_condition[[subtype_col]][nbr_idx],
      source_pt    = pt_i,
      target_pt    = df_condition[[pt_col]][nbr_idx],
      weight       = w,
      stringsAsFactors = FALSE
    )
  }
  
  edges <- bind_rows(edge_list)
  
  flux_mat <- matrix(
    0,
    nrow = length(state_levels),
    ncol = length(state_levels),
    dimnames = list(state_levels, state_levels)
  )
  
  if (nrow(edges) == 0) {
    return(list(edge_table = edges, flux_matrix = flux_mat))
  }
  
  source_counts <- df_condition %>%
    count(.data[[subtype_col]], name = "n_source_cells") %>%
    rename(source_state = 1)
  
  edge_sum <- edges %>%
    group_by(source_state, target_state) %>%
    summarize(weight_sum = sum(weight), .groups = "drop") %>%
    left_join(source_counts, by = "source_state") %>%
    mutate(flux = weight_sum / n_source_cells)
  
  for (ii in seq_len(nrow(edge_sum))) {
    flux_mat[edge_sum$source_state[ii], edge_sum$target_state[ii]] <- edge_sum$flux[ii]
  }
  
  list(edge_table = edges, flux_matrix = flux_mat)
}

# ============================================================
# 2. LOAD OBJECT
# ============================================================

load(obj_path)

if (!exists("s_objs")) {
  stop("Object 's_objs' was not found inside integrated_compartment_cancer_processed.RData")
}

obj <- s_objs

stopifnot(condition_col %in% colnames(obj@meta.data))
stopifnot(sample_col %in% colnames(obj@meta.data))
stopifnot(celltype_col %in% colnames(obj@meta.data))
stopifnot(umap_reduction %in% names(obj@reductions))

msg("Loaded object:")
print(obj)

msg("Condition table:")
print(table(obj@meta.data[[condition_col]], useNA = "ifany"))

msg("Cell-type table:")
print(table(obj@meta.data[[celltype_col]], useNA = "ifany"))

# ============================================================
# 3. ADD REGIMEN AND SUBSET TREATED SAMPLES
# ============================================================

meta0 <- obj@meta.data

meta0$chemo_regimen <- regimen_map[as.character(meta0[[sample_col]])]

obj$chemo_regimen <- meta0$chemo_regimen

# Make sure it was added
stopifnot("chemo_regimen" %in% colnames(obj@meta.data))

msg("Chemo regimen table before subsetting:")
print(table(obj@meta.data[[sample_col]], obj@meta.data$chemo_regimen, useNA = "ifany"))

cells_keep <- rownames(obj@meta.data)[
  obj@meta.data[[condition_col]] == treated_condition &
    !is.na(obj@meta.data$chemo_regimen) &
    obj@meta.data[[celltype_col]] %in% celltype_order
]

msg("Cells kept for treated-regimen trajectory: %s", length(cells_keep))

if (length(cells_keep) == 0) {
  stop("No cells passed the treated-regimen filter. Check condition values, sample IDs, and cell_type labels.")
}

obj_sub <- obj[, cells_keep]

obj_sub$chemo_regimen <- factor(
  obj_sub@meta.data$chemo_regimen,
  levels = regimen_levels
)

obj_sub$cell_type <- factor(
  obj_sub@meta.data[[celltype_col]],
  levels = celltype_order
)

msg("Treated sample x regimen table:")
print(table(obj_sub@meta.data[[sample_col]], obj_sub$chemo_regimen, useNA = "ifany"))

msg("Regimen x cell-type table:")
print(table(obj_sub$chemo_regimen, obj_sub$cell_type, useNA = "ifany"))

# ============================================================
# 4. STRATIFIED DOWNSAMPLING BY REGIMEN AND CELL TYPE
# ============================================================

meta <- obj_sub@meta.data
meta$cell <- rownames(meta)

sample_cells <- meta %>%
  group_by(chemo_regimen, cell_type) %>%
  group_modify(~{
    n_take <- min(nrow(.x), n_per_celltype_per_regimen)
    .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
  }) %>%
  ungroup()

obj_small <- obj_sub[, sample_cells$cell]

obj_small$chemo_regimen <- factor(
  obj_small$chemo_regimen,
  levels = regimen_levels
)

obj_small$cell_type <- factor(
  obj_small$cell_type,
  levels = celltype_order
)

msg("Downsampled regimen x cell-type table:")
print(table(obj_small$chemo_regimen, obj_small$cell_type, useNA = "ifany"))

# ============================================================
# 5. BUILD SHARED TREATED MONOCLE TRAJECTORY
# ============================================================

cds <- build_cds_with_seurat_umap(
  obj_small = obj_small,
  assay_counts = assay_counts,
  umap_reduction = umap_reduction,
  gene_min_cells = gene_min_cells,
  cluster_k = cluster_k,
  cluster_resolution = cluster_resolution,
  minimal_branch_len = minimal_branch_len
)

root_cells <- rownames(colData(cds))[
  as.character(colData(cds)[[celltype_col]]) == root_celltype
]

msg("Root cells found for %s: %s", root_celltype, length(root_cells))

if (length(root_cells) == 0) {
  stop(sprintf("No root cells found for root_celltype = '%s'", root_celltype))
}

if (length(root_cells) > 1000) {
  root_cells <- sample(root_cells, 1000)
}

cds <- order_cells(cds, root_cells = root_cells)

pt <- pseudotime(cds)
pt_scaled <- scale01(pt)
names(pt_scaled) <- names(pt)

obj_small$monocle_pseudotime <- NA_real_
obj_small$monocle_pseudotime[names(pt_scaled)] <- as.numeric(pt_scaled)

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.data.frame(closest_vertex)
closest_vertex$cell <- rownames(closest_vertex)
colnames(closest_vertex)[1] <- "graph_node"

# ============================================================
# 6. BUILD ANALYSIS DATAFRAME
# ============================================================

umap_mat <- Embeddings(obj_small, umap_reduction)
umap_df <- as.data.frame(umap_mat[, 1:2, drop = FALSE])
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)

meta_small <- obj_small@meta.data
meta_small$cell <- rownames(meta_small)

df <- umap_df %>%
  left_join(
    meta_small %>%
      transmute(
        cell = rownames(meta_small),
        cell_type = .data[[celltype_col]],
        chemo_regimen = chemo_regimen,
        sample_id = .data[[sample_col]]
      ),
    by = "cell"
  ) %>%
  left_join(
    tibble(cell = names(pt_scaled), pseudotime = as.numeric(pt_scaled)),
    by = "cell"
  ) %>%
  left_join(closest_vertex, by = "cell") %>%
  mutate(
    chemo_regimen = factor(chemo_regimen, levels = regimen_levels),
    cell_type = factor(cell_type, levels = celltype_order)
  )

write.csv(
  df,
  file.path(outdir, "Wilms_treated_regimen_cell_level_monocle_metadata.csv"),
  row.names = FALSE
)

# ============================================================
# 7. BASIC PLOTS
# ============================================================

p_regimen <- ggplot(df, aes(UMAP_1, UMAP_2, color = chemo_regimen)) +
  geom_point(size = 0.18, alpha = 0.75) +
  scale_color_manual(values = regimen_colors) +
  theme_classic() +
  labs(
    title = "Treated Wilms trajectory by chemotherapy regimen",
    color = "Regimen"
  )

save_gg(file.path(outdir, "UMAP_treated_regimen.pdf"), p_regimen, 8, 6.5)

p_celltype <- ggplot(df, aes(UMAP_1, UMAP_2, color = cell_type)) +
  geom_point(size = 0.18, alpha = 0.75) +
  theme_classic() +
  labs(
    title = "Treated Wilms shared trajectory by cell type",
    color = "Cell type"
  )

save_gg(file.path(outdir, "UMAP_treated_celltype.pdf"), p_celltype, 9, 7)

p_pt <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 0.18, alpha = 0.85, na.rm = TRUE) +
  scale_color_gradientn(
    colours = c("#2C0E7A", "#5B2A86", "#9A3C8E", "#D85C5C", "#F39B3D", "#F9E721"),
    limits = c(0, 1),
    oob = squish
  ) +
  theme_classic() +
  labs(
    title = "Shared treated Monocle3 pseudotime",
    color = "Pseudotime"
  )

save_gg(file.path(outdir, "UMAP_treated_pseudotime.pdf"), p_pt, 8, 6.5)

p_split_pt <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 0.18, alpha = 0.85, na.rm = TRUE) +
  facet_wrap(~chemo_regimen, nrow = 1) +
  scale_color_gradientn(
    colours = c("#2C0E7A", "#5B2A86", "#9A3C8E", "#D85C5C", "#F39B3D", "#F9E721"),
    limits = c(0, 1),
    oob = squish
  ) +
  theme_classic() +
  labs(
    title = "Shared treated pseudotime split by regimen",
    color = "Pseudotime"
  )

save_gg(file.path(outdir, "UMAP_treated_pseudotime_split_regimen.pdf"), p_split_pt, 13, 5.5)

p_box <- ggplot(df, aes(cell_type, pseudotime, fill = chemo_regimen)) +
  geom_boxplot(
    outlier.size = 0.15,
    width = 0.7,
    position = position_dodge(width = 0.8),
    na.rm = TRUE
  ) +
  scale_fill_manual(values = regimen_colors) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Treated Wilms pseudotime by cell type and regimen",
    x = NULL,
    y = "Scaled Monocle3 pseudotime",
    fill = "Regimen"
  )

save_gg(file.path(outdir, "Pseudotime_boxplot_celltype_regimen.pdf"), p_box, 15, 6)

# ============================================================
# 8. REGIMEN-SPECIFIC TRANSITION FLUX
# ============================================================

df_ee4a  <- df %>% filter(chemo_regimen == "EE-4A")
df_dd4a  <- df %>% filter(chemo_regimen == "DD-4A")
df_other <- df %>% filter(chemo_regimen == "Other")

flux_ee4a <- compute_directed_transition_flux(
  df_condition = df_ee4a,
  subtype_col = "cell_type",
  pt_col = "pseudotime",
  k_neighbors = k_neighbors,
  min_delta = pseudotime_forward_min_delta,
  max_neighbors_used = max_neighbors_used,
  state_levels = celltype_order
)

flux_dd4a <- compute_directed_transition_flux(
  df_condition = df_dd4a,
  subtype_col = "cell_type",
  pt_col = "pseudotime",
  k_neighbors = k_neighbors,
  min_delta = pseudotime_forward_min_delta,
  max_neighbors_used = max_neighbors_used,
  state_levels = celltype_order
)

flux_other <- compute_directed_transition_flux(
  df_condition = df_other,
  subtype_col = "cell_type",
  pt_col = "pseudotime",
  k_neighbors = k_neighbors,
  min_delta = pseudotime_forward_min_delta,
  max_neighbors_used = max_neighbors_used,
  state_levels = celltype_order
)

mat_ee4a  <- flux_ee4a$flux_matrix
mat_dd4a  <- flux_dd4a$flux_matrix
mat_other <- flux_other$flux_matrix

write.csv(mat_ee4a,  file.path(outdir, "Transition_flux_matrix_EE4A.csv"))
write.csv(mat_dd4a,  file.path(outdir, "Transition_flux_matrix_DD4A.csv"))
write.csv(mat_other, file.path(outdir, "Transition_flux_matrix_Other.csv"))

write.csv(
  matrix_to_long(mat_ee4a, "flux_EE4A"),
  file.path(outdir, "Transition_flux_matrix_EE4A_long.csv"),
  row.names = FALSE
)

write.csv(
  matrix_to_long(mat_dd4a, "flux_DD4A"),
  file.path(outdir, "Transition_flux_matrix_DD4A_long.csv"),
  row.names = FALSE
)

write.csv(
  matrix_to_long(mat_other, "flux_Other"),
  file.path(outdir, "Transition_flux_matrix_Other_long.csv"),
  row.names = FALSE
)

# Delta matrices
mat_EE4A_minus_DD4A  <- mat_ee4a - mat_dd4a
mat_EE4A_minus_Other <- mat_ee4a - mat_other
mat_DD4A_minus_Other <- mat_dd4a - mat_other

write.csv(mat_EE4A_minus_DD4A,  file.path(outdir, "Transition_flux_matrix_EE4A_minus_DD4A.csv"))
write.csv(mat_EE4A_minus_Other, file.path(outdir, "Transition_flux_matrix_EE4A_minus_Other.csv"))
write.csv(mat_DD4A_minus_Other, file.path(outdir, "Transition_flux_matrix_DD4A_minus_Other.csv"))

# Heatmaps
pdf(file.path(outdir, "Heatmap_transition_flux_EE4A.pdf"), width = 10, height = 9)
pheatmap(mat_ee4a, cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "grey85", main = "Directed transition flux: EE-4A")
dev.off()

pdf(file.path(outdir, "Heatmap_transition_flux_DD4A.pdf"), width = 10, height = 9)
pheatmap(mat_dd4a, cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "grey85", main = "Directed transition flux: DD-4A")
dev.off()

pdf(file.path(outdir, "Heatmap_transition_flux_Other.pdf"), width = 10, height = 9)
pheatmap(mat_other, cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "grey85", main = "Directed transition flux: Other")
dev.off()

pdf(file.path(outdir, "Heatmap_transition_flux_EE4A_minus_DD4A.pdf"), width = 10, height = 9)
pheatmap(mat_EE4A_minus_DD4A, cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "grey85", main = "Transition flux delta: EE-4A - DD-4A")
dev.off()

# ============================================================
# 9. EDGE TABLES
# ============================================================

edge_summary_ee4a <- flux_ee4a$edge_table %>%
  group_by(source_state, target_state) %>%
  summarize(weight_EE4A = sum(weight), .groups = "drop")

edge_summary_dd4a <- flux_dd4a$edge_table %>%
  group_by(source_state, target_state) %>%
  summarize(weight_DD4A = sum(weight), .groups = "drop")

edge_summary_other <- flux_other$edge_table %>%
  group_by(source_state, target_state) %>%
  summarize(weight_Other = sum(weight), .groups = "drop")

edge_compare <- full_join(edge_summary_ee4a, edge_summary_dd4a,
                          by = c("source_state", "target_state")) %>%
  full_join(edge_summary_other, by = c("source_state", "target_state")) %>%
  mutate(
    weight_EE4A = replace_na(weight_EE4A, 0),
    weight_DD4A = replace_na(weight_DD4A, 0),
    weight_Other = replace_na(weight_Other, 0),
    delta_EE4A_minus_DD4A = weight_EE4A - weight_DD4A,
    delta_EE4A_minus_Other = weight_EE4A - weight_Other,
    delta_DD4A_minus_Other = weight_DD4A - weight_Other,
    max_regimen = case_when(
      weight_EE4A >= weight_DD4A & weight_EE4A >= weight_Other ~ "EE-4A",
      weight_DD4A >= weight_EE4A & weight_DD4A >= weight_Other ~ "DD-4A",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(desc(abs(delta_EE4A_minus_DD4A)))

write.csv(
  edge_compare,
  file.path(outdir, "Edge_level_transition_flux_regimen_comparison.csv"),
  row.names = FALSE
)

top_edges <- edge_compare %>%
  filter(source_state != target_state) %>%
  arrange(desc(pmax(weight_EE4A, weight_DD4A, weight_Other))) %>%
  slice_head(n = 50)

write.csv(
  top_edges,
  file.path(outdir, "Top50_regimen_transition_edges.csv"),
  row.names = FALSE
)

# ============================================================
# 10. SHARED TREATED TRAJECTORY WITH REGIMEN-SPECIFIC ARROWS
# ============================================================

centroids <- df %>%
  filter(!is.na(cell_type), !is.na(UMAP_1), !is.na(UMAP_2), !is.na(pseudotime)) %>%
  group_by(cell_type) %>%
  summarize(
    x = median(UMAP_1),
    y = median(UMAP_2),
    n_cells = n(),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cell_type = as.character(cell_type))

edge_ee4a <- matrix_to_long(mat_ee4a, "flux") %>%
  filter(source != target, flux >= min_flux_to_draw) %>%
  mutate(chemo_regimen = "EE-4A")

edge_dd4a <- matrix_to_long(mat_dd4a, "flux") %>%
  filter(source != target, flux >= min_flux_to_draw) %>%
  mutate(chemo_regimen = "DD-4A")

edge_other <- matrix_to_long(mat_other, "flux") %>%
  filter(source != target, flux >= min_flux_to_draw) %>%
  mutate(chemo_regimen = "Other")

edges_overlay <- bind_rows(edge_ee4a, edge_dd4a, edge_other) %>%
  left_join(centroids %>% rename(source = cell_type, x = x, y = y), by = "source") %>%
  left_join(centroids %>% rename(target = cell_type, xend = x, yend = y), by = "target") %>%
  filter(!is.na(x), !is.na(y), !is.na(xend), !is.na(yend)) %>%
  group_by(chemo_regimen) %>%
  mutate(edge_width = scales::rescale(flux, to = c(0.7, 3.5))) %>%
  ungroup()

write.csv(
  edges_overlay,
  file.path(outdir, "Edges_used_for_regimen_arrow_overlay.csv"),
  row.names = FALSE
)

write.csv(
  centroids,
  file.path(outdir, "Celltype_centroids_for_regimen_arrow_overlay.csv"),
  row.names = FALSE
)

backbone_df <- centroids %>%
  arrange(median_pseudotime)

backbone_edges <- backbone_df %>%
  mutate(
    xend = lead(x),
    yend = lead(y),
    target = lead(cell_type)
  ) %>%
  rename(source = cell_type) %>%
  filter(!is.na(xend), !is.na(yend))

p_regimen_arrows <- ggplot() +
  geom_curve(
    data = backbone_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey40",
    linewidth = 0.7,
    linetype = "dotted",
    curvature = 0.1,
    alpha = 0.8
  ) +
  geom_curve(
    data = edges_overlay,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      linewidth = edge_width,
      color = chemo_regimen
    ),
    curvature = 0.18,
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    alpha = 0.9
  ) +
  geom_point(
    data = centroids,
    aes(x, y, size = n_cells),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 1.05
  ) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x, y, label = cell_type),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2
  ) +
  scale_color_manual(values = regimen_colors) +
  scale_linewidth_identity() +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  labs(
    title = "Treated Wilms trajectory by chemotherapy regimen",
    subtitle = "Dotted grey line = shared treated trajectory; red = EE-4A; green = DD-4A; brown = Other",
    color = "Regimen",
    size = "Cells"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right"
  )

save_gg(
  file.path(outdir, "Treated_shared_trajectory_EE4A_DD4A_Other_arrows.pdf"),
  p_regimen_arrows,
  width = 13,
  height = 10
)

# Version with grey background cells
p_regimen_arrows_background <- ggplot() +
  geom_point(
    data = df,
    aes(UMAP_1, UMAP_2),
    color = "grey85",
    size = 0.12,
    alpha = 0.25
  ) +
  geom_curve(
    data = backbone_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey40",
    linewidth = 0.7,
    linetype = "dotted",
    curvature = 0.1,
    alpha = 0.8
  ) +
  geom_curve(
    data = edges_overlay,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      linewidth = edge_width,
      color = chemo_regimen
    ),
    curvature = 0.18,
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    alpha = 0.9
  ) +
  geom_point(
    data = centroids,
    aes(x, y, size = n_cells),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 1.05
  ) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(x, y, label = cell_type),
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.2
  ) +
  scale_color_manual(values = regimen_colors) +
  scale_linewidth_identity() +
  scale_size_continuous(range = c(3, 10)) +
  theme_classic() +
  labs(
    title = "Treated Wilms trajectory by chemotherapy regimen",
    subtitle = "Dotted grey line = shared treated trajectory; red = EE-4A; green = DD-4A; brown = Other",
    color = "Regimen",
    size = "Cells"
  )

save_gg(
  file.path(outdir, "Treated_shared_trajectory_EE4A_DD4A_Other_arrows_with_background.pdf"),
  p_regimen_arrows_background,
  width = 13,
  height = 10
)

# ============================================================
# 11. NETWORK VERSION
# ============================================================

nodes_network <- centroids %>%
  transmute(
    name = cell_type,
    n_cells = n_cells,
    median_pseudotime = median_pseudotime
  )

edges_network <- edges_overlay %>%
  transmute(
    from = source,
    to = target,
    chemo_regimen = chemo_regimen,
    weight = flux
  )

g <- graph_from_data_frame(edges_network, vertices = nodes_network, directed = TRUE)

p_network <- ggraph(g, layout = "stress") +
  geom_edge_link(
    aes(width = weight, color = chemo_regimen),
    arrow = arrow(length = unit(3.5, "mm"), type = "closed"),
    end_cap = circle(5, "mm"),
    start_cap = circle(5, "mm"),
    alpha = 0.85
  ) +
  geom_node_point(aes(size = n_cells), shape = 21, fill = "white", color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.1) +
  scale_edge_color_manual(values = regimen_colors) +
  scale_edge_width(range = c(0.4, 3.2)) +
  scale_size_continuous(range = c(3, 10)) +
  theme_void() +
  labs(
    title = "Regimen-specific treated Wilms transition network",
    subtitle = "Red = EE-4A; green = DD-4A; brown = Other; edge width = transition flux"
  )

save_gg(
  file.path(outdir, "Network_regimen_transition_arrows.pdf"),
  p_network,
  width = 12,
  height = 9
)

# ============================================================
# 12. MONOCLE NATIVE PLOTS
# ============================================================

pdf(file.path(outdir, "Monocle_plot_cells_pseudotime_treated_regimen.pdf"), width = 8, height = 7)
print(
  plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE
  )
)
dev.off()

pdf(file.path(outdir, "Monocle_plot_cells_celltype_treated_regimen.pdf"), width = 9, height = 7)
print(
  plot_cells(
    cds,
    color_cells_by = celltype_col,
    label_cell_groups = TRUE,
    label_groups_by_cluster = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE
  )
)
dev.off()

pdf(file.path(outdir, "Monocle_plot_cells_chemo_regimen.pdf"), width = 8, height = 7)
print(
  plot_cells(
    cds,
    color_cells_by = "chemo_regimen",
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE
  )
)
dev.off()

# ============================================================
# 13. SAVE OBJECTS
# ============================================================

saveRDS(
  obj_small,
  file.path(outdir, "Wilms_treated_regimen_monocle_seurat_subset.rds")
)

saveRDS(
  cds,
  file.path(outdir, "Wilms_treated_regimen_monocle_cds.rds")
)

saveRDS(
  df,
  file.path(outdir, "Wilms_treated_regimen_monocle_dataframe.rds")
)

msg("Analysis completed.")
msg("Outputs written to: %s", normalizePath(outdir))
