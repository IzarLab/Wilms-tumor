# ============================================================
# Figure 3 untreated malignant subcompartment trajectory
# ============================================================
#
# GitHub script name : 01_analyze_untreated_subcompartment_monocle_trajectory.R
# Original file      : Fig.3.monocle.wilms.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Runs Monocle3 trajectory analysis in untreated malignant cells using blastema as root and compares blastema, epithelial, and stromal states.
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

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(scales)
  library(patchwork)
  library(pheatmap)
  library(igraph)
  library(ggraph)
})

set.seed(1234)

# ============================================================
# 0. SETTINGS
# ============================================================

obj_path <- "/home/ubuntu/wilms/integrated_compartment_cancer_processed.RData"
outdir   <- "./Wilms_untreated_subcompartment_monocle_transition"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

condition_col <- "condition"
state_col     <- "subcompartment"
sample_col    <- "orig.ident"

condition_keep <- "untreated"

states_keep <- c("blastema", "stroma", "epithelium")
state_order <- c("blastema", "epithelium", "stroma")

root_state <- "blastema"

assay_counts <- "RNA"
umap_reduction <- "umap"

n_per_state <- 4000

gene_min_cells <- 10
cluster_k <- 20
cluster_resolution <- 1e-3
minimal_branch_len <- 10

k_neighbors <- 20
pseudotime_forward_min_delta <- 0.01
max_neighbors_used <- 8

n_pt_bins <- 8

# ============================================================
# 1. HELPERS
# ============================================================

msg <- function(...) cat(sprintf(...), "\n")

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(NA_real_, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

get_counts_safe <- function(obj, assay = "RNA") {
  
  assays_available <- names(obj@assays)
  
  if (is.null(assays_available) || length(assays_available) == 0) {
    stop("No assays found in obj@assays.")
  }
  
  if (!assay %in% assays_available) {
    stop(sprintf(
      "Assay '%s' not found. Available assays: %s",
      assay,
      paste(assays_available, collapse = ", ")
    ))
  }
  
  mat <- NULL
  
  # Seurat v5 layer syntax
  mat <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  
  # Seurat v4 slot syntax
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
  ggsave(filename, plot = plot, width = width, height = height,
         dpi = 300, device = cairo_pdf)
}

get_knn_indices <- function(mat, k = 20) {
  d <- as.matrix(dist(mat))
  diag(d) <- Inf
  idx <- t(apply(d, 1, function(z) order(z)[seq_len(min(k, length(z)))]))
  idx
}

matrix_to_long <- function(mat, value_name = "value") {
  as.data.frame(as.table(mat), stringsAsFactors = FALSE) %>%
    rename(source = Var1, target = Var2, !!value_name := Freq)
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
    subtype_col = "state",
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
stopifnot(state_col %in% colnames(obj@meta.data))
stopifnot(sample_col %in% colnames(obj@meta.data))
stopifnot(umap_reduction %in% names(obj@reductions))

msg("Original object:")
print(obj)

msg("Condition table:")
print(table(obj@meta.data[[condition_col]], useNA = "ifany"))

msg("Subcompartment table:")
print(table(obj@meta.data[[state_col]], useNA = "ifany"))

# ============================================================
# 3. KEEP ONLY UNTREATED WILMS SUBCOMPARTMENTS
# ============================================================

obj_sub <- subset(
  obj,
  subset = condition == condition_keep & subcompartment %in% states_keep
)

obj_sub[[state_col]] <- factor(obj_sub@meta.data[[state_col]], levels = state_order)

msg("Subset object:")
print(obj_sub)

msg("Untreated subcompartment counts:")
print(table(obj_sub@meta.data[[state_col]], useNA = "ifany"))

# Downsample for Monocle speed / balanced trajectory learning
meta <- obj_sub@meta.data
meta$cell <- rownames(meta)

sample_cells <- meta %>%
  group_by(.data[[state_col]]) %>%
  group_modify(~{
    n_take <- min(nrow(.x), n_per_state)
    .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
  }) %>%
  ungroup()

obj_small <- obj_sub[, sample_cells$cell]

msg("Downsampled counts:")
print(table(obj_small@meta.data[[state_col]], useNA = "ifany"))

# ============================================================
# 4. BUILD MONOCLE3 CDS USING RNA COUNTS + EXISTING SEURAT UMAP
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

root_cells <- rownames(colData(cds))[as.character(colData(cds)[[state_col]]) == root_state]

if (length(root_cells) == 0) {
  stop(sprintf("No root cells found for root_state = '%s'", root_state))
}

if (length(root_cells) > 500) {
  root_cells <- sample(root_cells, 500)
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
# 5. ANALYSIS DATAFRAME
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
        state = .data[[state_col]],
        condition = .data[[condition_col]],
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
    state = factor(state, levels = state_order)
  )

write.csv(df, file.path(outdir, "Wilms_untreated_cell_level_monocle_metadata.csv"),
          row.names = FALSE)

# ============================================================
# 6. BASIC UMAP / PSEUDOTIME PLOTS
# ============================================================

p_state <- ggplot(df, aes(UMAP_1, UMAP_2, color = state)) +
  geom_point(size = 0.25, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Wilms untreated subcompartments",
    color = "Subcompartment"
  )

save_gg(file.path(outdir, "UMAP_untreated_subcompartment.pdf"), p_state, 7, 6)

p_pt <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
  geom_point(size = 0.25, alpha = 0.85, na.rm = TRUE) +
  scale_color_gradientn(
    colours = c("#2C0E7A", "#5B2A86", "#9A3C8E", "#D85C5C", "#F39B3D", "#F9E721"),
    limits = c(0, 1),
    oob = squish
  ) +
  theme_classic() +
  labs(
    title = "Wilms untreated Monocle3 pseudotime",
    color = "Pseudotime"
  )

save_gg(file.path(outdir, "UMAP_untreated_monocle_pseudotime.pdf"), p_pt, 7, 6)

p_box <- ggplot(df, aes(state, pseudotime, fill = state)) +
  geom_boxplot(outlier.size = 0.2, width = 0.7, na.rm = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(
    title = "Pseudotime by Wilms subcompartment",
    x = NULL,
    y = "Scaled Monocle3 pseudotime"
  ) +
  guides(fill = "none")

save_gg(file.path(outdir, "Pseudotime_boxplot_by_subcompartment.pdf"), p_box, 6.5, 5)

p_density <- ggplot(df, aes(pseudotime, fill = state, color = state)) +
  geom_density(alpha = 0.18, na.rm = TRUE) +
  theme_classic() +
  labs(
    title = "Pseudotime density by Wilms subcompartment",
    x = "Scaled Monocle3 pseudotime",
    y = "Density"
  )

save_gg(file.path(outdir, "Pseudotime_density_by_subcompartment.pdf"), p_density, 7, 5)

# ============================================================
# 7. PSEUDOTIME BIN COMPOSITION
# ============================================================

df_bins <- df %>%
  filter(!is.na(pseudotime)) %>%
  mutate(
    pt_bin = cut(
      pseudotime,
      breaks = seq(0, 1, length.out = n_pt_bins + 1),
      include.lowest = TRUE,
      labels = paste0("Bin", seq_len(n_pt_bins))
    )
  )

bin_comp <- df_bins %>%
  count(pt_bin, state, name = "n_cells") %>%
  group_by(pt_bin) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(bin_comp, file.path(outdir, "Pseudotime_bin_subcompartment_composition.csv"),
          row.names = FALSE)

p_bins <- ggplot(bin_comp, aes(pt_bin, frac, fill = state)) +
  geom_col(position = "fill") +
  theme_classic() +
  labs(
    title = "Wilms subcompartment composition across pseudotime",
    x = "Pseudotime bin",
    y = "Relative composition",
    fill = "Subcompartment"
  )

save_gg(file.path(outdir, "Pseudotime_bin_subcompartment_composition.pdf"), p_bins, 8, 5)

# ============================================================
# 8. DIRECTED STATE TRANSITION FLUX
# ============================================================

flux_res <- compute_directed_transition_flux(
  df_condition = df,
  subtype_col = "state",
  pt_col = "pseudotime",
  k_neighbors = k_neighbors,
  min_delta = pseudotime_forward_min_delta,
  max_neighbors_used = max_neighbors_used,
  state_levels = state_order
)

mat_flux <- flux_res$flux_matrix

write.csv(mat_flux, file.path(outdir, "Transition_flux_matrix_untreated.csv"))
write.csv(matrix_to_long(mat_flux, "flux"),
          file.path(outdir, "Transition_flux_matrix_untreated_long.csv"),
          row.names = FALSE)

edge_summary <- flux_res$edge_table %>%
  group_by(source_state, target_state) %>%
  summarize(weight = sum(weight), .groups = "drop") %>%
  arrange(desc(weight))

write.csv(edge_summary, file.path(outdir, "Edge_level_transition_summary_untreated.csv"),
          row.names = FALSE)

pdf(file.path(outdir, "Heatmap_transition_flux_untreated.pdf"), width = 6.5, height = 5.8)
pheatmap(
  mat_flux,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey85",
  main = "Directed transition flux: untreated Wilms"
)
dev.off()

# ============================================================
# 9. ADDITIONAL FIGURE: STATE TRANSITIONS WITH ARROWS
# ============================================================

edge_df <- matrix_to_long(mat_flux, "flux") %>%
  filter(source != target, flux > 0) %>%
  arrange(desc(flux))

write.csv(edge_df, file.path(outdir, "Transition_edges_for_arrow_plot.csv"),
          row.names = FALSE)

centroids <- df %>%
  filter(!is.na(state), !is.na(UMAP_1), !is.na(UMAP_2)) %>%
  group_by(state) %>%
  summarize(
    x = median(UMAP_1),
    y = median(UMAP_2),
    n_cells = n(),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(state = as.character(state))

edge_plot_df <- edge_df %>%
  left_join(centroids %>% rename(source = state, x = x, y = y), by = "source") %>%
  left_join(centroids %>% rename(target = state, xend = x, yend = y), by = "target") %>%
  mutate(
    edge_width = scales::rescale(flux, to = c(0.6, 3.5)),
    label = sprintf("%.3f", flux)
  )

p_arrow_umap <- ggplot() +
  geom_point(
    data = df,
    aes(UMAP_1, UMAP_2, color = state),
    size = 0.18,
    alpha = 0.22
  ) +
  geom_curve(
    data = edge_plot_df,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = edge_width),
    curvature = 0.18,
    arrow = arrow(length = unit(0.18, "inches"), type = "closed"),
    color = "black",
    alpha = 0.85
  ) +
  geom_point(
    data = centroids,
    aes(x, y, size = n_cells),
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 1.1
  ) +
  geom_text(
    data = centroids,
    aes(x, y, label = state),
    size = 4,
    fontface = "bold",
    vjust = -1
  ) +
  scale_linewidth_identity() +
  scale_size_continuous(range = c(5, 12)) +
  theme_classic() +
  labs(
    title = "Untreated Wilms state-transition arrows on UMAP",
    subtitle = "Arrows point from earlier to later pseudotime neighbors; width reflects transition flux",
    color = "Subcompartment",
    size = "Cells"
  )

save_gg(file.path(outdir, "UMAP_state_transition_arrows_untreated.pdf"),
        p_arrow_umap, 8, 7)

# Network-only version
nodes <- centroids %>%
  transmute(
    name = state,
    n_cells = n_cells,
    median_pseudotime = median_pseudotime
  )

edges_network <- edge_df %>%
  rename(from = source, to = target, weight = flux)

g <- graph_from_data_frame(edges_network, vertices = nodes, directed = TRUE)

p_network <- ggraph(g, layout = "circle") +
  geom_edge_link(
    aes(width = weight),
    arrow = arrow(length = unit(4, "mm"), type = "closed"),
    end_cap = circle(6, "mm"),
    start_cap = circle(6, "mm"),
    alpha = 0.85
  ) +
  geom_node_point(aes(size = n_cells), shape = 21, fill = "white", color = "black") +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.5, 3.5)) +
  scale_size_continuous(range = c(6, 13)) +
  theme_void() +
  labs(
    title = "Untreated Wilms subcompartment transition network",
    subtitle = "Directed edges are inferred from local forward pseudotime neighborhoods"
  )

save_gg(file.path(outdir, "Network_state_transition_arrows_untreated.pdf"),
        p_network, 7, 6)

# ============================================================
# 10. MONOCLE NATIVE PLOTS
# ============================================================

pdf(file.path(outdir, "Monocle_plot_cells_pseudotime_untreated.pdf"), width = 8, height = 7)
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

pdf(file.path(outdir, "Monocle_plot_cells_subcompartment_untreated.pdf"), width = 8, height = 7)
print(
  plot_cells(
    cds,
    color_cells_by = state_col,
    label_cell_groups = TRUE,
    label_groups_by_cluster = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE
  )
)
dev.off()

# ============================================================
# 11. SAVE OBJECTS
# ============================================================

saveRDS(obj_small, file.path(outdir, "Wilms_untreated_subcompartment_monocle_seurat_subset.rds"))
saveRDS(cds, file.path(outdir, "Wilms_untreated_subcompartment_monocle_cds.rds"))
saveRDS(df, file.path(outdir, "Wilms_untreated_subcompartment_monocle_dataframe.rds"))

msg("Analysis completed.")
msg("Outputs written to: %s", normalizePath(outdir))
