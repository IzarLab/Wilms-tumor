# ============================================================
# Joint treated versus untreated malignant trajectory analysis
# ============================================================
#
# GitHub script name : 01_analyze_joint_treated_vs_untreated_monocle_trajectory.R
# Original file      : Fig.3.monocle.wilms.joint.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Runs joint Monocle3 trajectory analysis comparing untreated and neoadjuvant-treated malignant cell states.
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

outdir <- "./Wilms_joint_treated_vs_untreated_celltype_monocle"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

condition_col <- "condition"
celltype_col  <- "cell_type"
sample_col    <- "orig.ident"

assay_counts <- "RNA"
umap_reduction <- "umap"

condition_levels <- c("untreated", "neoadjuvant")

root_celltype <- "blastema 1"

# Cell-type order based on your table
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

# Downsampling
n_per_celltype_per_condition <- 1200

# Monocle parameters
gene_min_cells <- 10
cluster_k <- 20
cluster_resolution <- 1e-3
minimal_branch_len <- 10

# Transition parameters
k_neighbors <- 20
pseudotime_forward_min_delta <- 0.01
max_neighbors_used <- 8
min_flux_to_draw <- 0.005

# Plot colors
untreated_arrow_col <- "#2166AC"
  treated_arrow_col   <- "#B2182B"
    trajectory_col      <- "grey35"
      
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
    
    extract_monocle_graph_edges <- function(cds) {
      g <- principal_graph(cds)[["UMAP"]]
      coords <- t(as.matrix(igraph::layout_as_tree(g)))
      
      graph_layout <- principal_graph_aux(cds)[["UMAP"]]$dp_mst
      
      if (!is.null(graph_layout)) {
        node_coords <- as.data.frame(t(graph_layout))
        node_coords$node <- rownames(node_coords)
        colnames(node_coords)[1:2] <- c("x", "y")
      } else {
        node_coords <- data.frame()
      }
      
      edge_df <- igraph::as_data_frame(g, what = "edges")
      
      if (nrow(node_coords) > 0) {
        edge_df <- edge_df %>%
          left_join(node_coords %>% rename(from = node, x = x, y = y), by = "from") %>%
          left_join(node_coords %>% rename(to = node, xend = x, yend = y), by = "to")
      }
      
      edge_df
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
    stopifnot(celltype_col %in% colnames(obj@meta.data))
    stopifnot(sample_col %in% colnames(obj@meta.data))
    stopifnot(umap_reduction %in% names(obj@reductions))
    
    msg("Loaded object:")
    print(obj)
    
    msg("Condition table:")
    print(table(obj@meta.data[[condition_col]], useNA = "ifany"))
    
    msg("Cell-type table:")
    print(table(obj@meta.data[[celltype_col]], useNA = "ifany"))
    
    # ============================================================
    # 3. SUBSET TREATED + UNTREATED CELL TYPES
    # ============================================================
    
    obj_sub <- subset(
      obj,
      subset = condition %in% condition_levels & cell_type %in% celltype_order
    )
    
    obj_sub[[condition_col]] <- factor(
      obj_sub@meta.data[[condition_col]],
      levels = condition_levels
    )
    
    obj_sub[[celltype_col]] <- factor(
      obj_sub@meta.data[[celltype_col]],
      levels = celltype_order
    )
    
    msg("Subset condition x cell_type table:")
    print(table(obj_sub@meta.data[[condition_col]], obj_sub@meta.data[[celltype_col]]))
    
    # ============================================================
    # 4. STRATIFIED DOWNSAMPLING
    # ============================================================
    
    meta <- obj_sub@meta.data
    meta$cell <- rownames(meta)
    
    sample_cells <- meta %>%
      group_by(.data[[condition_col]], .data[[celltype_col]]) %>%
      group_modify(~{
        n_take <- min(nrow(.x), n_per_celltype_per_condition)
        .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
      }) %>%
      ungroup()
    
    obj_small <- obj_sub[, sample_cells$cell]
    
    obj_small[[condition_col]] <- factor(
      obj_small@meta.data[[condition_col]],
      levels = condition_levels
    )
    
    obj_small[[celltype_col]] <- factor(
      obj_small@meta.data[[celltype_col]],
      levels = celltype_order
    )
    
    msg("Downsampled condition x cell_type table:")
    print(table(obj_small@meta.data[[condition_col]], obj_small@meta.data[[celltype_col]]))
    
    # ============================================================
    # 5. BUILD JOINT SHARED MONOCLE TRAJECTORY
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
    
    # Closest graph vertex
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.data.frame(closest_vertex)
    closest_vertex$cell <- rownames(closest_vertex)
    colnames(closest_vertex)[1] <- "graph_node"
    
    # ============================================================
    # 6. ANALYSIS DATAFRAME
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
        condition = factor(condition, levels = condition_levels),
        cell_type = factor(cell_type, levels = celltype_order)
      )
    
    write.csv(
      df,
      file.path(outdir, "Wilms_joint_cell_level_monocle_metadata.csv"),
      row.names = FALSE
    )
    
    # ============================================================
    # 7. BASIC SHARED TRAJECTORY PLOTS
    # ============================================================
    
    p_celltype <- ggplot(df, aes(UMAP_1, UMAP_2, color = cell_type)) +
      geom_point(size = 0.18, alpha = 0.75) +
      theme_classic() +
      labs(
        title = "Joint Wilms trajectory by cell type",
        color = "Cell type"
      )
    
    save_gg(file.path(outdir, "UMAP_joint_celltype.pdf"), p_celltype, 9, 7)
    
    p_condition <- ggplot(df, aes(UMAP_1, UMAP_2, color = condition)) +
      geom_point(size = 0.18, alpha = 0.7) +
      scale_color_manual(
        values = c(
          untreated = untreated_arrow_col,
          neoadjuvant = treated_arrow_col
        )
      ) +
      theme_classic() +
      labs(
        title = "Joint Wilms trajectory by treatment condition",
        color = "Condition"
      )
    
    save_gg(file.path(outdir, "UMAP_joint_condition.pdf"), p_condition, 8, 6.5)
    
    p_pt <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
      geom_point(size = 0.18, alpha = 0.85, na.rm = TRUE) +
      scale_color_gradientn(
        colours = c("#2C0E7A", "#5B2A86", "#9A3C8E", "#D85C5C", "#F39B3D", "#F9E721"),
        limits = c(0, 1),
        oob = squish
      ) +
      theme_classic() +
      labs(
        title = "Joint Wilms Monocle3 pseudotime",
        color = "Pseudotime"
      )
    
    save_gg(file.path(outdir, "UMAP_joint_pseudotime.pdf"), p_pt, 8, 6.5)
    
    p_split_pt <- ggplot(df, aes(UMAP_1, UMAP_2, color = pseudotime)) +
      geom_point(size = 0.18, alpha = 0.85, na.rm = TRUE) +
      facet_wrap(~condition, nrow = 1) +
      scale_color_gradientn(
        colours = c("#2C0E7A", "#5B2A86", "#9A3C8E", "#D85C5C", "#F39B3D", "#F9E721"),
        limits = c(0, 1),
        oob = squish
      ) +
      theme_classic() +
      labs(
        title = "Joint pseudotime split by condition",
        color = "Pseudotime"
      )
    
    save_gg(file.path(outdir, "UMAP_joint_pseudotime_split_condition.pdf"), p_split_pt, 12, 5.5)
    
    # ============================================================
    # 8. PSEUDOTIME SUMMARY
    # ============================================================
    
    pt_stats_celltype <- df %>%
      filter(!is.na(pseudotime), !is.na(cell_type), !is.na(condition)) %>%
      group_by(cell_type) %>%
      summarize(
        n_untreated = sum(condition == "untreated"),
        n_treated = sum(condition == "neoadjuvant"),
        median_untreated = median(pseudotime[condition == "untreated"], na.rm = TRUE),
        median_treated = median(pseudotime[condition == "neoadjuvant"], na.rm = TRUE),
        delta_treated_minus_untreated = median_treated - median_untreated,
        p_value = tryCatch(
          wilcox.test(
            pseudotime[condition == "untreated"],
            pseudotime[condition == "neoadjuvant"]
          )$p.value,
          error = function(e) NA_real_
        ),
        .groups = "drop"
      ) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(p_adj, desc(abs(delta_treated_minus_untreated)))
    
    write.csv(
      pt_stats_celltype,
      file.path(outdir, "Pseudotime_stats_treated_vs_untreated_by_celltype.csv"),
      row.names = FALSE
    )
    
    p_box <- ggplot(df, aes(cell_type, pseudotime, fill = condition)) +
      geom_boxplot(
        outlier.size = 0.15,
        width = 0.7,
        position = position_dodge(width = 0.8),
        na.rm = TRUE
      ) +
      scale_fill_manual(
        values = c(
          untreated = untreated_arrow_col,
          neoadjuvant = treated_arrow_col
        )
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Joint pseudotime by cell type and condition",
        x = NULL,
        y = "Scaled Monocle3 pseudotime",
        fill = "Condition"
      )
    
    save_gg(file.path(outdir, "Pseudotime_boxplot_celltype_condition.pdf"), p_box, 15, 6)
    
    # ============================================================
    # 9. CONDITION-SPECIFIC TRANSITION FLUX ON SHARED TRAJECTORY
    # ============================================================
    
    df_untreated <- df %>% filter(condition == "untreated")
    df_treated   <- df %>% filter(condition == "neoadjuvant")
    
    flux_untreated <- compute_directed_transition_flux(
      df_condition = df_untreated,
      subtype_col = "cell_type",
      pt_col = "pseudotime",
      k_neighbors = k_neighbors,
      min_delta = pseudotime_forward_min_delta,
      max_neighbors_used = max_neighbors_used,
      state_levels = celltype_order
    )
    
    flux_treated <- compute_directed_transition_flux(
      df_condition = df_treated,
      subtype_col = "cell_type",
      pt_col = "pseudotime",
      k_neighbors = k_neighbors,
      min_delta = pseudotime_forward_min_delta,
      max_neighbors_used = max_neighbors_used,
      state_levels = celltype_order
    )
    
    mat_untreated <- flux_untreated$flux_matrix
    mat_treated   <- flux_treated$flux_matrix
    mat_delta     <- mat_treated - mat_untreated
    mat_log2fc    <- log2((mat_treated + 1e-6) / (mat_untreated + 1e-6))
    
    write.csv(mat_untreated, file.path(outdir, "Transition_flux_matrix_untreated.csv"))
    write.csv(mat_treated, file.path(outdir, "Transition_flux_matrix_treated.csv"))
    write.csv(mat_delta, file.path(outdir, "Transition_flux_matrix_treated_minus_untreated.csv"))
    write.csv(mat_log2fc, file.path(outdir, "Transition_flux_matrix_log2FC_treated_vs_untreated.csv"))
    
    write.csv(
      matrix_to_long(mat_untreated, "flux_untreated"),
      file.path(outdir, "Transition_flux_matrix_untreated_long.csv"),
      row.names = FALSE
    )
    
    write.csv(
      matrix_to_long(mat_treated, "flux_treated"),
      file.path(outdir, "Transition_flux_matrix_treated_long.csv"),
      row.names = FALSE
    )
    
    write.csv(
      matrix_to_long(mat_delta, "delta_flux"),
      file.path(outdir, "Transition_flux_matrix_delta_long.csv"),
      row.names = FALSE
    )
    
    # Heatmaps
    pdf(file.path(outdir, "Heatmap_transition_flux_untreated.pdf"), width = 10, height = 9)
    pheatmap(
      mat_untreated,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "grey85",
      main = "Directed transition flux: untreated"
    )
    dev.off()
    
    pdf(file.path(outdir, "Heatmap_transition_flux_treated.pdf"), width = 10, height = 9)
    pheatmap(
      mat_treated,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "grey85",
      main = "Directed transition flux: treated/neoadjuvant"
    )
    dev.off()
    
    pdf(file.path(outdir, "Heatmap_transition_flux_treated_minus_untreated.pdf"), width = 10, height = 9)
    pheatmap(
      mat_delta,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "grey85",
      main = "Transition flux delta: treated - untreated"
    )
    dev.off()
    
    # ============================================================
    # 10. PREPARE ARROW OVERLAY DATA
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
    
    edge_untreated <- matrix_to_long(mat_untreated, "flux") %>%
      filter(source != target, flux >= min_flux_to_draw) %>%
      mutate(condition = "untreated")
    
    edge_treated <- matrix_to_long(mat_treated, "flux") %>%
      filter(source != target, flux >= min_flux_to_draw) %>%
      mutate(condition = "neoadjuvant")
    
    edges_overlay <- bind_rows(edge_untreated, edge_treated) %>%
      left_join(centroids %>% rename(source = cell_type, x = x, y = y), by = "source") %>%
      left_join(centroids %>% rename(target = cell_type, xend = x, yend = y), by = "target") %>%
      filter(!is.na(x), !is.na(y), !is.na(xend), !is.na(yend)) %>%
      group_by(condition) %>%
      mutate(edge_width = scales::rescale(flux, to = c(0.4, 2.8))) %>%
      ungroup()
    
    write.csv(
      edges_overlay,
      file.path(outdir, "Edges_used_for_joint_arrow_overlay.csv"),
      row.names = FALSE
    )
    
    write.csv(
      centroids,
      file.path(outdir, "Celltype_centroids_for_arrow_overlay.csv"),
      row.names = FALSE
    )
    
    # ============================================================
    # 11. SHARED TRAJECTORY DOTTED + CONDITION ARROWS
    # ============================================================
    
    # Approximate shared backbone using centroid order by pseudotime
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
    
    p_joint_arrows <- ggplot() +
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
        color = trajectory_col,
        linewidth = 0.55,
        linetype = "dotted",
        curvature = 0.1,
        alpha = 0.75
      ) +
      geom_curve(
        data = edges_overlay,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend,
          linewidth = edge_width,
          color = condition
        ),
        curvature = 0.18,
        arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
        alpha = 0.82
      ) +
      geom_point(
        data = centroids,
        aes(x, y, size = n_cells),
        shape = 21,
        fill = "white",
        color = "black",
        stroke = 1.05
      ) +
      geom_text_repel(
        data = centroids,
        aes(x, y, label = cell_type),
        size = 3.0,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2
      ) +
      scale_color_manual(
        values = c(
          untreated = untreated_arrow_col,
          neoadjuvant = treated_arrow_col
        ),
        labels = c(
          untreated = "Untreated",
          neoadjuvant = "Treated/neoadjuvant"
        )
      ) +
      scale_linewidth_identity() +
      scale_size_continuous(range = c(3, 10)) +
      theme_classic() +
      labs(
        title = "Integrated Wilms trajectory with condition-specific state transitions",
        subtitle = "Dotted line = shared joint trajectory; blue arrows = untreated; red arrows = treated/neoadjuvant",
        color = "Condition",
        size = "Cells"
      )
    
    save_gg(
      file.path(outdir, "Joint_shared_trajectory_treated_vs_untreated_arrows.pdf"),
      p_joint_arrows,
      width = 13,
      height = 10
    )
    
    # Same figure without grey background cells
    p_joint_arrows_clean <- ggplot() +
      geom_curve(
        data = backbone_edges,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = trajectory_col,
        linewidth = 0.65,
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
          color = condition
        ),
        curvature = 0.18,
        arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
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
      geom_text_repel(
        data = centroids,
        aes(x, y, label = cell_type),
        size = 3.2,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2
      ) +
      scale_color_manual(
        values = c(
          untreated = untreated_arrow_col,
          neoadjuvant = treated_arrow_col
        ),
        labels = c(
          untreated = "Untreated",
          neoadjuvant = "Treated/neoadjuvant"
        )
      ) +
      scale_linewidth_identity() +
      scale_size_continuous(range = c(3, 10)) +
      theme_void() +
      labs(
        title = "Integrated Wilms state-transition map",
        subtitle = "Dotted line = shared joint trajectory; blue arrows = untreated; red arrows = treated/neoadjuvant",
        color = "Condition",
        size = "Cells"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        legend.position = "right"
      )
    
    save_gg(
      file.path(outdir, "Joint_shared_trajectory_treated_vs_untreated_arrows_clean.pdf"),
      p_joint_arrows_clean,
      width = 13,
      height = 10
    )
    
    # ============================================================
    # 12. NETWORK VERSION
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
        condition = condition,
        weight = flux
      )
    
    g <- graph_from_data_frame(edges_network, vertices = nodes_network, directed = TRUE)
    
    p_network <- ggraph(g, layout = "stress") +
      geom_edge_link(
        aes(width = weight, color = condition),
        arrow = arrow(length = unit(3.5, "mm"), type = "closed"),
        end_cap = circle(5, "mm"),
        start_cap = circle(5, "mm"),
        alpha = 0.82
      ) +
      geom_node_point(aes(size = n_cells), shape = 21, fill = "white", color = "black") +
      geom_node_text(aes(label = name), repel = TRUE, size = 3.1) +
      scale_edge_color_manual(
        values = c(
          untreated = untreated_arrow_col,
          neoadjuvant = treated_arrow_col
        ),
        labels = c(
          untreated = "Untreated",
          neoadjuvant = "Treated/neoadjuvant"
        )
      ) +
      scale_edge_width(range = c(0.4, 2.8)) +
      scale_size_continuous(range = c(3, 10)) +
      theme_void() +
      labs(
        title = "Condition-specific Wilms transition network",
        subtitle = "Blue = untreated; red = treated/neoadjuvant; edge width = transition flux"
      )
    
    save_gg(
      file.path(outdir, "Network_treated_vs_untreated_transition_arrows.pdf"),
      p_network,
      width = 12,
      height = 9
    )
    
    # ============================================================
    # 13. EDGE COMPARISON TABLE
    # ============================================================
    
    edge_summary_untreated <- flux_untreated$edge_table %>%
      group_by(source_state, target_state) %>%
      summarize(weight_untreated = sum(weight), .groups = "drop")
    
    edge_summary_treated <- flux_treated$edge_table %>%
      group_by(source_state, target_state) %>%
      summarize(weight_treated = sum(weight), .groups = "drop")
    
    edge_compare <- full_join(
      edge_summary_untreated,
      edge_summary_treated,
      by = c("source_state", "target_state")
    ) %>%
      mutate(
        weight_untreated = replace_na(weight_untreated, 0),
        weight_treated = replace_na(weight_treated, 0),
        delta_treated_minus_untreated = weight_treated - weight_untreated,
        log2FC_treated_vs_untreated = log2((weight_treated + 1e-6) / (weight_untreated + 1e-6))
      ) %>%
      arrange(desc(abs(delta_treated_minus_untreated)))
    
    write.csv(
      edge_compare,
      file.path(outdir, "Edge_level_transition_rewiring_treated_vs_untreated.csv"),
      row.names = FALSE
    )
    
    top_edges <- edge_compare %>%
      filter(source_state != target_state) %>%
      slice_head(n = 30)
    
    write.csv(
      top_edges,
      file.path(outdir, "Top30_rewired_transitions_treated_vs_untreated.csv"),
      row.names = FALSE
    )
    
    # ============================================================
    # 14. MONOCLE NATIVE PLOTS
    # ============================================================
    
    pdf(file.path(outdir, "Monocle_plot_cells_pseudotime_joint.pdf"), width = 8, height = 7)
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
    
    pdf(file.path(outdir, "Monocle_plot_cells_celltype_joint.pdf"), width = 9, height = 7)
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
    
    pdf(file.path(outdir, "Monocle_plot_cells_condition_joint.pdf"), width = 8, height = 7)
    print(
      plot_cells(
        cds,
        color_cells_by = condition_col,
        label_cell_groups = FALSE,
        label_groups_by_cluster = FALSE,
        label_branch_points = TRUE,
        label_leaves = TRUE
      )
    )
    dev.off()
    
    # ============================================================
    # 15. SAVE OBJECTS
    # ============================================================
    
    saveRDS(
      obj_small,
      file.path(outdir, "Wilms_joint_treated_untreated_monocle_seurat_subset.rds")
    )
    
    saveRDS(
      cds,
      file.path(outdir, "Wilms_joint_treated_untreated_monocle_cds.rds")
    )
    
    saveRDS(
      df,
      file.path(outdir, "Wilms_joint_treated_untreated_monocle_dataframe.rds")
    )
    
    msg("Analysis completed.")
    msg("Outputs written to: %s", normalizePath(outdir))