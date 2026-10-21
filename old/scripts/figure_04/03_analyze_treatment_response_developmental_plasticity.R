# ============================================================
# Treatment-response developmental plasticity analysis
# ============================================================
#
# GitHub script name : 03_analyze_treatment_response_developmental_plasticity.R
# Original file      : Treatment_regimes.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Analyzes developmental plasticity, composition, response groups, and treatment-regimen-level remodeling using current Wilms TME object annotations.
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

# # ============================================================
# # Wilms tumor developmental plasticity vs treatment response
# # Uses current object structure:
# # Phase, subcompartment, compartment_name, condition, cell_type, cell_type_2
# # ============================================================
# 
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(dplyr)
#   library(tidyr)
#   library(tibble)
#   library(ggplot2)
#   library(ggpubr)
#   library(patchwork)
#   library(scales)
#   library(forcats)
#   library(stringr)
# })
# 
# # ============================================================
# # Input / output
# # ============================================================
# 
# infile <- "~/Documents/Izar_Group/wilms/TME.RData"
# outdir <- "~/Documents/Izar_Group/wilms/Plasticity_Response_Analysis_CurrentObject"
# 
# dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# 
# load(infile)
# obj <- tme_
# 
# message("Loaded object with ", ncol(obj), " cells")
# 
# # ============================================================
# # Inspect metadata
# # ============================================================
# 
# print(colnames(obj@meta.data))
# print(table(obj@meta.data$condition))
# print(table(obj@meta.data$subcompartment))
# print(table(obj@meta.data$cell_type_2))
# 
# # ============================================================
# # Sample metadata from clinical table
# # ============================================================
# 
# sample_meta <- tibble::tribble(
#   ~study_id, ~sample, ~tumor_analyzed_code, ~neoadjuvant_code, ~response_code,
#   ~cell_type_code, ~histology_code, ~stage, ~recurrence_code,
#   "EC015", "W612", 1, "N/A", "N/A", "1",   1, 2, "1",
#   "EG010", "W002", 1, "3",   "1",   "2",   2, 4, "1",
#   "EG015", "W340", 1, "N/A", "N/A", "2",   1, 4, "1",
#   "EC013", "W011", 1, "N/A", "N/A", "2",   1, 3, "2",
#   "EG013", "W050", 1, "5*",  "2",   "1",   1, 3, "1",
#   "EG047", "W010", 1, "1",   "1",   "6",   1, 2, "2",
#   "EG048", "W944", 1, "2",   "1",   "1,5", 1, 4, "2",
#   "EC011", "W132", 1, "N/A", "N/A", "1",   1, 4, "2",
#   "EG044", "W009", 1, "1",   "3",   "4",   1, 5, "2",
#   "EG043", "W184", 1, "2",   "1",   "1",   2, 3, "2",
#   "EG042", "W495", 1, "N/A", "N/A", "3",   1, 2, "2",
#   "EG041", "W008", 1, "N/A", "N/A", "1",   1, 3, "2",
#   "EC007", "W007", 1, "N/A", "N/A", "1",   1, 2, "2",
#   "EG049", "W180", 2, "2",   "2",   "2",   1, 3, "1",
#   "EG037", "W370", 1, "2",   "1",   "1",   2, 5, "2",
#   "EC016", "W006", 2, "2",   "UNK", "1",   1, 3, "1",
#   "EG031", "W005", 1, "N/A", "N/A", "2",   1, 2, "2",
#   "EG027", "W004", 1, "N/A", "N/A", "2",   1, 1, "1",
#   "EG026", "W241", 1, "N/A", "N/A", "1",   2, 3, "UNK",
#   "EG023", "W040", 1, "N/A", "N/A", "1",   1, 2, "1",
#   "EG022", "W003", 1, "1",   "1",   "1",   1, 4, "1",
#   "EG004", "W001", 2, "UNK", "UNK", "4",   1, 1, "1",
#   "EG002", "W441", 1, "1",   "2",   "2",   1, 4, "2",
#   "EG051", "W012", 1, "N/A", "N/A", "2",   1, 3, "2"
# )
# 
# sample_meta <- sample_meta %>%
#   mutate(
#     tumor_analyzed = case_when(
#       tumor_analyzed_code == 1 ~ "Primary",
#       tumor_analyzed_code == 2 ~ "Recurrence",
#       TRUE ~ as.character(tumor_analyzed_code)
#     ),
#     treatment = case_when(
#       neoadjuvant_code == "1" ~ "EE-4A",
#       neoadjuvant_code == "2" ~ "DD-4A",
#       neoadjuvant_code == "3" ~ "I",
#       neoadjuvant_code == "4" ~ "M",
#       neoadjuvant_code %in% c("5", "5*") ~ "V+D",
#       neoadjuvant_code == "N/A" ~ "Untreated",
#       neoadjuvant_code == "UNK" ~ "UNK",
#       TRUE ~ neoadjuvant_code
#     ),
#     response = case_when(
#       response_code == "1" ~ "Dec",
#       response_code == "2" ~ "Stable",
#       response_code == "3" ~ "Inc",
#       response_code == "N/A" ~ "N/A",
#       response_code == "UNK" ~ "UNK",
#       TRUE ~ response_code
#     ),
#     response_group = case_when(
#       response == "Dec" ~ "Responder",
#       response %in% c("Stable", "Inc") ~ "Non-responder",
#       TRUE ~ NA_character_
#     ),
#     histology = case_when(
#       histology_code == 1 ~ "Favorable",
#       histology_code == 2 ~ "Anaplastic",
#       TRUE ~ as.character(histology_code)
#     ),
#     recurrence = case_when(
#       recurrence_code == "1" ~ "Yes",
#       recurrence_code == "2" ~ "No",
#       recurrence_code == "UNK" ~ "UNK",
#       TRUE ~ recurrence_code
#     ),
#     stage = as.character(stage)
#   )
# 
# write.csv(sample_meta, file.path(outdir, "sample_metadata_cleaned.csv"), row.names = FALSE)
# 
# # ============================================================
# # Harmonize sample IDs in object
# # ============================================================
# 
# obj$sample <- as.character(obj$orig.ident)
# 
# obj$sample <- case_when(
#   obj$sample == "W1"  ~ "W001",
#   obj$sample == "W2"  ~ "W002",
#   obj$sample == "W3"  ~ "W003",
#   obj$sample == "W4"  ~ "W004",
#   obj$sample == "W5"  ~ "W005",
#   obj$sample == "W6"  ~ "W006",
#   obj$sample == "W7"  ~ "W007",
#   obj$sample == "W8"  ~ "W008",
#   obj$sample == "W9"  ~ "W009",
#   obj$sample == "W10" ~ "W010",
#   obj$sample == "W11" ~ "W011",
#   obj$sample == "W12" ~ "W012",
#   obj$sample == "W40" ~ "W040",
#   obj$sample == "W50" ~ "W050",
#   TRUE ~ obj$sample
# )
# 
# # ============================================================
# # Add sample metadata to object
# # ============================================================
# 
# meta_add <- sample_meta %>%
#   select(
#     sample, study_id, tumor_analyzed, treatment, response,
#     response_group, histology, stage, recurrence
#   )
# 
# obj@meta.data <- obj@meta.data %>%
#   rownames_to_column("cell_id") %>%
#   left_join(meta_add, by = "sample") %>%
#   column_to_rownames("cell_id")
# 
# message("Samples after metadata merge:")
# print(table(obj$sample, useNA = "ifany"))
# 
# message("Response group distribution:")
# print(table(obj$response_group, useNA = "ifany"))
# 
# # ============================================================
# # Use refined cell annotation
# # ============================================================
# 
# state_col <- "cell_type_2"
# 
# obj$state <- as.character(obj@meta.data[[state_col]])
# 
# # ============================================================
# # Define developmental categories
# # ============================================================
# 
# blastema_states <- c(
#   "blastema 1",
#   "blastema 2",
#   "blastema 3",
#   "S blastema 1",
#   "mitotic blastema 1",
#   "mitotic blastema 2",
#   "G2M blastema 1",
#   "G2M blastema 2"
# )
# 
# epithelial_states <- c(
#   "tubules",
#   "mitotic tubules",
#   "G2M tubules",
#   "S tubules",
#   "podocyte-like"
# )
# 
# stromal_myogenic_states <- c(
#   "smooth myocyte-like",
#   "smooth myocyte-like precursors",
#   "nascent smooth myocyte-like",
#   "mitotic smooth myocyte-like",
#   "mitotic smooth myocyte-like precursors",
#   "G2M smooth myocyte-like",
#   "G2M smooth myocyte-like precursors",
#   "S smooth myocyte-like",
#   "PAX3+ myogenic precursors",
#   "myoblast-like",
#   "striated myocyte-like"
# )
# 
# obj$developmental_group <- case_when(
#   obj$state %in% blastema_states ~ "Blastema",
#   obj$state %in% epithelial_states ~ "Epithelial",
#   obj$state %in% stromal_myogenic_states ~ "Stromal/Myogenic",
#   obj$subcompartment == "stroma" ~ "Stromal/Myogenic",
#   obj$subcompartment == "epithelium" ~ "Epithelial",
#   obj$subcompartment == "blastema" ~ "Blastema",
#   TRUE ~ "Other"
# )
# 
# obj$developmental_group <- factor(
#   obj$developmental_group,
#   levels = c("Blastema", "Epithelial", "Stromal/Myogenic", "Other")
# )
# 
# message("Developmental group counts:")
# print(table(obj$developmental_group, useNA = "ifany"))
# 
# # ============================================================
# # Restrict to malignant cells only
# # ============================================================
# 
# obj.mal <- subset(
#   obj,
#   subset = compartment_name == "cancer"
# )
# 
# message("Malignant cells:")
# print(table(obj.mal$developmental_group))
# 
# # ============================================================
# # Treated samples with known response
# # ============================================================
# 
# obj.tr <- subset(
#   obj.mal,
#   subset = response_group %in% c("Responder", "Non-responder")
# )
# 
# message("Treated known-response malignant cells:")
# print(table(obj.tr$response_group))
# print(table(obj.tr$response_group, obj.tr$developmental_group))
# 
# # ============================================================
# # Sample-level developmental composition
# # ============================================================
# 
# prop_df <- obj.tr@meta.data %>%
#   count(
#     sample, response_group, response, treatment, histology,
#     stage, recurrence, developmental_group,
#     name = "n"
#   ) %>%
#   group_by(sample) %>%
#   mutate(
#     total_malignant = sum(n),
#     prop = n / total_malignant
#   ) %>%
#   ungroup()
# 
# write.csv(
#   prop_df,
#   file.path(outdir, "sample_level_developmental_group_proportions.csv"),
#   row.names = FALSE
# )
# 
# # ============================================================
# # State-level sample composition
# # ============================================================
# 
# state_prop_df <- obj.tr@meta.data %>%
#   count(
#     sample, response_group, response, treatment, histology,
#     stage, recurrence, state,
#     name = "n"
#   ) %>%
#   group_by(sample) %>%
#   mutate(
#     total_malignant = sum(n),
#     prop = n / total_malignant
#   ) %>%
#   ungroup()
# 
# write.csv(
#   state_prop_df,
#   file.path(outdir, "sample_level_cell_state_proportions.csv"),
#   row.names = FALSE
# )
# 
# # ============================================================
# # Plasticity / resistance metrics per sample
# # ============================================================
# 
# metric_df <- prop_df %>%
#   select(
#     sample, response_group, response, treatment, histology,
#     stage, recurrence, developmental_group, prop
#   ) %>%
#   pivot_wider(
#     names_from = developmental_group,
#     values_from = prop,
#     values_fill = 0
#   ) %>%
#   mutate(
#     blastema_persistence_score = Blastema,
#     epithelial_differentiation_score = Epithelial,
#     stromal_myogenic_differentiation_score = `Stromal/Myogenic`,
#     total_differentiation_score = Epithelial + `Stromal/Myogenic`,
#     plasticity_ratio = total_differentiation_score / (blastema_persistence_score + 1e-6),
#     resistance_index = blastema_persistence_score / (total_differentiation_score + 1e-6)
#   )
# 
# # Developmental entropy
# entropy_df <- prop_df %>%
#   group_by(sample, response_group, response, treatment, histology, stage, recurrence) %>%
#   summarise(
#     developmental_entropy = -sum(prop * log(prop + 1e-9)),
#     .groups = "drop"
#   )
# 
# metric_df <- metric_df %>%
#   left_join(
#     entropy_df,
#     by = c("sample", "response_group", "response", "treatment",
#            "histology", "stage", "recurrence")
#   )
# 
# # Proliferation among blastema
# blastema_prolif_df <- obj.tr@meta.data %>%
#   filter(developmental_group == "Blastema") %>%
#   count(sample, Phase, name = "n") %>%
#   group_by(sample) %>%
#   mutate(prop = n / sum(n)) %>%
#   ungroup() %>%
#   filter(Phase %in% c("S", "G2M")) %>%
#   group_by(sample) %>%
#   summarise(
#     cycling_blastema_fraction = sum(prop),
#     .groups = "drop"
#   )
# 
# metric_df <- metric_df %>%
#   left_join(blastema_prolif_df, by = "sample") %>%
#   mutate(
#     cycling_blastema_fraction = ifelse(
#       is.na(cycling_blastema_fraction),
#       0,
#       cycling_blastema_fraction
#     )
#   )
# 
# write.csv(
#   metric_df,
#   file.path(outdir, "plasticity_response_metrics.csv"),
#   row.names = FALSE
# )
# 
# # ============================================================
# # Statistics: responder vs non-responder
# # ============================================================
# 
# test_vars <- c(
#   "blastema_persistence_score",
#   "total_differentiation_score",
#   "epithelial_differentiation_score",
#   "stromal_myogenic_differentiation_score",
#   "plasticity_ratio",
#   "resistance_index",
#   "developmental_entropy",
#   "cycling_blastema_fraction"
# )
# 
# stats_df <- lapply(test_vars, function(v) {
#   
#   x <- metric_df %>%
#     filter(
#       !is.na(.data[[v]]),
#       response_group %in% c("Responder", "Non-responder")
#     )
#   
#   if (length(unique(x$response_group)) < 2) {
#     return(tibble(metric = v, p_value = NA_real_))
#   }
#   
#   test <- wilcox.test(
#     as.formula(paste(v, "~ response_group")),
#     data = x,
#     exact = FALSE
#   )
#   
#   tibble(
#     metric = v,
#     p_value = test$p.value
#   )
#   
# }) %>%
#   bind_rows() %>%
#   mutate(
#     p_adj = p.adjust(p_value, method = "BH"),
#     significance = case_when(
#       is.na(p_adj) ~ NA_character_,
#       p_adj < 0.0001 ~ "****",
#       p_adj < 0.001  ~ "***",
#       p_adj < 0.01   ~ "**",
#       p_adj < 0.05   ~ "*",
#       TRUE ~ "ns"
#     )
#   )
# 
# write.csv(
#   stats_df,
#   file.path(outdir, "responder_vs_nonresponder_statistics.csv"),
#   row.names = FALSE
# )
# 
# print(stats_df)
# 
# # ============================================================
# # Plot colors
# # ============================================================
# 
# resp_cols <- c(
#   "Responder" = "#4DAF4A",
#   "Non-responder" = "#E41A1C"
# )
# 
# dev_cols <- c(
#   "Blastema" = "#E64B35",
#   "Epithelial" = "#4DBBD5",
#   "Stromal/Myogenic" = "#00A087",
#   "Other" = "grey80"
# )
# 
# # ============================================================
# # Plot 1: sample-wise malignant developmental composition
# # ============================================================
# 
# sample_order <- metric_df %>%
#   arrange(response_group, desc(blastema_persistence_score)) %>%
#   pull(sample)
# 
# prop_df$sample <- factor(prop_df$sample, levels = sample_order)
# 
# p_sample_stack <- ggplot(
#   prop_df,
#   aes(x = sample, y = prop * 100, fill = developmental_group)
# ) +
#   geom_col(color = "black", linewidth = 0.2) +
#   facet_grid(. ~ response_group, scales = "free_x", space = "free_x") +
#   scale_fill_manual(values = dev_cols, drop = FALSE) +
#   labs(
#     x = NULL,
#     y = "Malignant composition (%)",
#     fill = "Developmental group",
#     title = "Responders and non-responders differ in malignant developmental-state composition"
#   ) +
#   theme_classic(base_size = 13) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   file.path(outdir, "01_samplewise_developmental_composition_response.pdf"),
#   p_sample_stack,
#   width = 10,
#   height = 5
# )
# 
# # ============================================================
# # Helper for metric plots
# # ============================================================
# 
# plot_metric <- function(df, metric, ylab, title) {
#   pval <- stats_df %>%
#     filter(metric == !!metric) %>%
#     pull(significance)
#   
#   ggplot(
#     df,
#     aes(x = response_group, y = .data[[metric]], fill = response_group)
#   ) +
#     geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
#     geom_jitter(width = 0.12, size = 3, alpha = 0.9) +
#     scale_fill_manual(values = resp_cols) +
#     labs(
#       x = NULL,
#       y = ylab,
#       title = paste0(title, " (", pval, ")")
#     ) +
#     theme_classic(base_size = 13) +
#     theme(
#       legend.position = "none",
#       plot.title = element_text(face = "bold", size = 11)
#     )
# }
# 
# p1 <- plot_metric(
#   metric_df,
#   "blastema_persistence_score",
#   "Blastema fraction",
#   "Blastema persistence"
# )
# 
# p2 <- plot_metric(
#   metric_df,
#   "total_differentiation_score",
#   "Differentiated fraction",
#   "Total differentiation"
# )
# 
# p3 <- plot_metric(
#   metric_df,
#   "plasticity_ratio",
#   "Differentiated / Blastema",
#   "Plasticity ratio"
# )
# 
# p4 <- plot_metric(
#   metric_df,
#   "resistance_index",
#   "Blastema / Differentiated",
#   "Resistance index"
# )
# 
# p5 <- plot_metric(
#   metric_df,
#   "developmental_entropy",
#   "Developmental entropy",
#   "State diversity"
# )
# 
# p6 <- plot_metric(
#   metric_df,
#   "cycling_blastema_fraction",
#   "Cycling blastema fraction",
#   "Cycling blastema"
# )
# 
# combined_metrics <- (p1 | p2 | p3) / (p4 | p5 | p6)
# 
# ggsave(
#   file.path(outdir, "02_response_plasticity_metrics.pdf"),
#   combined_metrics,
#   width = 13,
#   height = 8
# )
# 
# # ============================================================
# # Plot 3: cell-state heatmap by response
# # ============================================================
# 
# state_order <- c(
#   "blastema 1",
#   "S blastema 1",
#   "mitotic blastema 1",
#   "blastema 2",
#   "mitotic blastema 2",
#   "blastema 3",
#   "tubules",
#   "S tubules",
#   "mitotic tubules",
#   "podocyte-like",
#   "PAX3+ myogenic precursors",
#   "myoblast-like",
#   "striated myocyte-like",
#   "smooth myocyte-like precursors",
#   "mitotic smooth myocyte-like precursors",
#   "nascent smooth myocyte-like",
#   "smooth myocyte-like",
#   "mitotic smooth myocyte-like"
# )
# 
# state_prop_df <- state_prop_df %>%
#   mutate(
#     sample = factor(sample, levels = sample_order),
#     state = factor(state, levels = rev(state_order))
#   )
# 
# p_state_heat <- ggplot(
#   state_prop_df,
#   aes(x = sample, y = state, fill = prop)
# ) +
#   geom_tile(color = "white", linewidth = 0.2) +
#   facet_grid(. ~ response_group, scales = "free_x", space = "free_x") +
#   scale_fill_gradient(
#     low = "white",
#     high = "black",
#     labels = percent_format()
#   ) +
#   labs(
#     x = NULL,
#     y = NULL,
#     fill = "Fraction",
#     title = "Fine malignant cell-state composition by treatment response"
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     axis.text.y = element_text(size = 9),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   file.path(outdir, "03_cell_state_composition_heatmap_response.pdf"),
#   p_state_heat,
#   width = 12,
#   height = 7
# )
# 
# # ============================================================
# # Plot 4: treatment-regimen stratified plasticity
# # ============================================================
# 
# metric_df$treatment <- factor(
#   metric_df$treatment,
#   levels = c("EE-4A", "DD-4A", "V+D", "I", "UNK")
# )
# 
# p_treatment <- ggplot(
#   metric_df,
#   aes(x = treatment, y = plasticity_ratio, fill = response_group)
# ) +
#   geom_boxplot(
#     outlier.shape = NA,
#     alpha = 0.7,
#     position = position_dodge(width = 0.8)
#   ) +
#   geom_jitter(
#     aes(color = response_group),
#     position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
#     size = 3
#   ) +
#   scale_fill_manual(values = resp_cols, na.value = "grey80") +
#   scale_color_manual(values = resp_cols, na.value = "grey80") +
#   labs(
#     x = "Treatment regimen",
#     y = "Plasticity ratio\n(Differentiated / Blastema)",
#     fill = "Response",
#     color = "Response",
#     title = "Plasticity ratio by treatment regimen and response"
#   ) +
#   theme_classic(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   file.path(outdir, "04_plasticity_ratio_by_treatment_response.pdf"),
#   p_treatment,
#   width = 8,
#   height = 5
# )
# 
# # ============================================================
# # Plot 5: response + histology
# # ============================================================
# 
# p_histology <- ggplot(
#   metric_df,
#   aes(x = histology, y = resistance_index, fill = response_group)
# ) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
#   geom_jitter(
#     aes(color = response_group),
#     position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
#     size = 3
#   ) +
#   scale_fill_manual(values = resp_cols, na.value = "grey80") +
#   scale_color_manual(values = resp_cols, na.value = "grey80") +
#   labs(
#     x = "Histology",
#     y = "Resistance index\n(Blastema / Differentiated)",
#     fill = "Response",
#     color = "Response",
#     title = "Blastema persistence and resistance by histology"
#   ) +
#   theme_classic(base_size = 13)
# 
# ggsave(
#   file.path(outdir, "05_resistance_index_by_histology_response.pdf"),
#   p_histology,
#   width = 7,
#   height = 5
# )
# 
# # ============================================================
# # Optional UMAPs
# # ============================================================
# 
# umap_reduction <- if ("umap" %in% names(obj@reductions)) {
#   "umap"
# } else if ("umap.recomputed" %in% names(obj@reductions)) {
#   "umap.recomputed"
# } else {
#   NA_character_
# }
# 
# if (!is.na(umap_reduction)) {
#   
#   p_umap_resp <- DimPlot(
#     obj.tr,
#     reduction = umap_reduction,
#     group.by = "response_group",
#     cols = resp_cols,
#     pt.size = 0.15,
#     raster = FALSE
#   ) +
#     theme_classic(base_size = 13) +
#     ggtitle("Malignant cells by treatment response")
#   
#   ggsave(
#     file.path(outdir, "06_UMAP_response_group.pdf"),
#     p_umap_resp,
#     width = 6,
#     height = 5
#   )
#   
#   p_umap_dev <- DimPlot(
#     obj.tr,
#     reduction = umap_reduction,
#     group.by = "developmental_group",
#     cols = dev_cols,
#     pt.size = 0.15,
#     raster = FALSE
#   ) +
#     theme_classic(base_size = 13) +
#     ggtitle("Malignant developmental groups")
#   
#   ggsave(
#     file.path(outdir, "07_UMAP_developmental_group.pdf"),
#     p_umap_dev,
#     width = 6,
#     height = 5
#   )
# }
# 
# # ============================================================
# # Save final object and summary
# # ============================================================
# 
# saveRDS(
#   obj,
#   file.path(outdir, "Wilms_TME_with_response_plasticity_metadata.rds")
# )
# 
# message("Done.")
# message("Outputs written to: ", outdir)

# ============================================================
# Wilms tumor developmental plasticity vs treatment response
# Uses current object structure:
# Phase, subcompartment, compartment_name, condition, cell_type, cell_type_2
# ============================================================

# ============================================================
# Wilms tumor developmental plasticity vs treatment response
# UPDATED response definition:
#   Include only treated samples with response_code 1/2/3
#   1 = Responder (decreased tumor burden)
#   2/3 = Non-responder (unchanged/increased tumor burden)
#
# Uses current object structure:
#   Phase, subcompartment, compartment_name, condition, cell_type, cell_type_2
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(forcats)
  library(stringr)
})

# ============================================================
# Input / output
# ============================================================

infile <- "~/Documents/Izar_Group/wilms/TME.RData"
outdir <- "~/Documents/Izar_Group/wilms/Plasticity_Response_Analysis_Response123"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

load(infile)
obj <- tme_

message("Loaded object with ", ncol(obj), " cells")

# ============================================================
# Sample metadata from clinical table
# ============================================================

sample_meta <- tibble::tribble(
  ~study_id, ~sample, ~tumor_analyzed_code, ~neoadjuvant_code, ~response_code,
  ~cell_type_code, ~histology_code, ~stage, ~recurrence_code,
  "EC015", "W612", 1, "N/A", "N/A", "1",   1, 2, "1",
  "EG010", "W002", 1, "3",   "1",   "2",   2, 4, "1",
  "EG015", "W340", 1, "N/A", "N/A", "2",   1, 4, "1",
  "EC013", "W011", 1, "N/A", "N/A", "2",   1, 3, "2",
  "EG013", "W050", 1, "5*",  "2",   "1",   1, 3, "1",
  "EG047", "W010", 1, "1",   "1",   "6",   1, 2, "2",
  "EG048", "W944", 1, "2",   "1",   "1,5", 1, 4, "2",
  "EC011", "W132", 1, "N/A", "N/A", "1",   1, 4, "2",
  "EG044", "W009", 1, "1",   "3",   "4",   1, 5, "2",
  "EG043", "W184", 1, "2",   "1",   "1",   2, 3, "2",
  "EG042", "W495", 1, "N/A", "N/A", "3",   1, 2, "2",
  "EG041", "W008", 1, "N/A", "N/A", "1",   1, 3, "2",
  "EC007", "W007", 1, "N/A", "N/A", "1",   1, 2, "2",
  "EG049", "W180", 2, "2",   "2",   "2",   1, 3, "1",
  "EG037", "W370", 1, "2",   "1",   "1",   2, 5, "2",
  "EC016", "W006", 2, "2",   "UNK", "1",   1, 3, "1",
  "EG031", "W005", 1, "N/A", "N/A", "2",   1, 2, "2",
  "EG027", "W004", 1, "N/A", "N/A", "2",   1, 1, "1",
  "EG026", "W241", 1, "N/A", "N/A", "1",   2, 3, "UNK",
  "EG023", "W040", 1, "N/A", "N/A", "1",   1, 2, "1",
  "EG022", "W003", 1, "1",   "1",   "1",   1, 4, "1",
  "EG004", "W001", 2, "UNK", "UNK", "4",   1, 1, "1",
  "EG002", "W441", 1, "1",   "2",   "2",   1, 4, "2",
  "EG051", "W012", 1, "N/A", "N/A", "2",   1, 3, "2"
)

sample_meta <- sample_meta %>%
  mutate(
    tumor_analyzed = case_when(
      tumor_analyzed_code == 1 ~ "Primary",
      tumor_analyzed_code == 2 ~ "Recurrence",
      TRUE ~ as.character(tumor_analyzed_code)
    ),
    
    treatment = case_when(
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
      response_code == "1" ~ "Decreased tumor burden",
      response_code == "2" ~ "Unchanged tumor burden",
      response_code == "3" ~ "Increased tumor burden",
      response_code == "N/A" ~ "N/A",
      response_code == "UNK" ~ "UNK",
      TRUE ~ response_code
    ),
    
    response_group = case_when(
      response_code == "1" ~ "Responder",
      response_code %in% c("2", "3") ~ "Non-responder",
      TRUE ~ NA_character_
    ),
    
    response_binary = case_when(
      response_code == "1" ~ 1,
      response_code %in% c("2", "3") ~ 0,
      TRUE ~ NA_real_
    ),
    
    histology = case_when(
      histology_code == 1 ~ "Favorable",
      histology_code == 2 ~ "Anaplastic",
      TRUE ~ as.character(histology_code)
    ),
    
    recurrence = case_when(
      recurrence_code == "1" ~ "Yes",
      recurrence_code == "2" ~ "No",
      recurrence_code == "UNK" ~ "UNK",
      TRUE ~ recurrence_code
    ),
    
    stage = as.character(stage)
  )

write.csv(
  sample_meta,
  file.path(outdir, "sample_metadata_cleaned_response123.csv"),
  row.names = FALSE
)

# ============================================================
# Harmonize sample IDs
# ============================================================

obj$sample <- as.character(obj$orig.ident)

obj$sample <- case_when(
  obj$sample == "W1"  ~ "W001",
  obj$sample == "W2"  ~ "W002",
  obj$sample == "W3"  ~ "W003",
  obj$sample == "W4"  ~ "W004",
  obj$sample == "W5"  ~ "W005",
  obj$sample == "W6"  ~ "W006",
  obj$sample == "W7"  ~ "W007",
  obj$sample == "W8"  ~ "W008",
  obj$sample == "W9"  ~ "W009",
  obj$sample == "W10" ~ "W010",
  obj$sample == "W11" ~ "W011",
  obj$sample == "W12" ~ "W012",
  obj$sample == "W40" ~ "W040",
  obj$sample == "W50" ~ "W050",
  TRUE ~ obj$sample
)

# ============================================================
# Add sample metadata
# ============================================================

meta_add <- sample_meta %>%
  select(
    sample,
    study_id,
    tumor_analyzed,
    treatment,
    response,
    response_code,
    response_group,
    response_binary,
    histology,
    stage,
    recurrence
  )

obj@meta.data <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(meta_add, by = "sample") %>%
  column_to_rownames("cell_id")

message("Response groups after merge:")
print(table(obj$response_group, useNA = "ifany"))

# ============================================================
# Annotation column
# ============================================================

state_col <- "cell_type_2"
obj$state <- as.character(obj@meta.data[[state_col]])

# ============================================================
# Define malignant developmental categories
# ============================================================

blastema_states <- c(
  "blastema 1",
  "blastema 2",
  "blastema 3",
  "S blastema 1",
  "mitotic blastema 1",
  "mitotic blastema 2",
  "G2M blastema 1",
  "G2M blastema 2"
)

epithelial_states <- c(
  "tubules",
  "mitotic tubules",
  "G2M tubules",
  "S tubules",
  "podocyte-like"
)

stromal_myogenic_states <- c(
  "smooth myocyte-like",
  "smooth myocyte-like precursors",
  "nascent smooth myocyte-like",
  "mitotic smooth myocyte-like",
  "mitotic smooth myocyte-like precursors",
  "G2M smooth myocyte-like",
  "G2M smooth myocyte-like precursors",
  "S smooth myocyte-like",
  "PAX3+ myogenic precursors",
  "myoblast-like",
  "striated myocyte-like"
)

obj$developmental_group <- case_when(
  obj$state %in% blastema_states ~ "Blastema",
  obj$state %in% epithelial_states ~ "Epithelial",
  obj$state %in% stromal_myogenic_states ~ "Stromal/Myogenic",
  obj$subcompartment == "stroma" ~ "Stromal/Myogenic",
  obj$subcompartment == "epithelium" ~ "Epithelial",
  obj$subcompartment == "blastema" ~ "Blastema",
  TRUE ~ "Other"
)

obj$developmental_group <- factor(
  obj$developmental_group,
  levels = c("Blastema", "Epithelial", "Stromal/Myogenic", "Other")
)

# ============================================================
# Keep malignant cells only
# ============================================================

obj.mal <- subset(
  obj,
  subset = compartment_name == "cancer"
)

message("Malignant developmental groups:")
print(table(obj.mal$developmental_group, useNA = "ifany"))

# ============================================================
# Keep treated samples with response 1/2/3 only
# ============================================================

obj.tr <- subset(
  obj.mal,
  subset = response_code %in% c("1", "2", "3")
)

message("Treated malignant cells with response 1/2/3:")
print(table(obj.tr$sample, obj.tr$response_group))
print(table(obj.tr$response_group, obj.tr$developmental_group))

included_samples <- obj.tr@meta.data %>%
  distinct(
    sample, response_code, response, response_group,
    treatment, histology, stage, recurrence
  ) %>%
  arrange(response_group, sample)

write.csv(
  included_samples,
  file.path(outdir, "included_samples_response123.csv"),
  row.names = FALSE
)

# ============================================================
# Sample-level developmental composition
# ============================================================

prop_df <- obj.tr@meta.data %>%
  count(
    sample,
    response_code,
    response,
    response_group,
    response_binary,
    treatment,
    histology,
    stage,
    recurrence,
    developmental_group,
    name = "n"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_malignant = sum(n),
    prop = n / total_malignant
  ) %>%
  ungroup()

write.csv(
  prop_df,
  file.path(outdir, "sample_level_developmental_group_proportions_response123.csv"),
  row.names = FALSE
)

# ============================================================
# Fine state-level composition
# ============================================================

state_prop_df <- obj.tr@meta.data %>%
  count(
    sample,
    response_code,
    response,
    response_group,
    response_binary,
    treatment,
    histology,
    stage,
    recurrence,
    state,
    name = "n"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_malignant = sum(n),
    prop = n / total_malignant
  ) %>%
  ungroup()

write.csv(
  state_prop_df,
  file.path(outdir, "sample_level_cell_state_proportions_response123.csv"),
  row.names = FALSE
)

# ============================================================
# Plasticity / resistance metrics
# ============================================================

metric_df <- prop_df %>%
  select(
    sample,
    response_code,
    response,
    response_group,
    response_binary,
    treatment,
    histology,
    stage,
    recurrence,
    developmental_group,
    prop
  ) %>%
  pivot_wider(
    names_from = developmental_group,
    values_from = prop,
    values_fill = 0
  ) %>%
  mutate(
    blastema_persistence_score = Blastema,
    epithelial_differentiation_score = Epithelial,
    stromal_myogenic_differentiation_score = `Stromal/Myogenic`,
    total_differentiation_score = Epithelial + `Stromal/Myogenic`,
    
    # Use log-ratio to avoid huge inflation when blastema is near zero
    plasticity_log_ratio = log2((total_differentiation_score + 0.01) /
                                  (blastema_persistence_score + 0.01)),
    
    resistance_log_ratio = log2((blastema_persistence_score + 0.01) /
                                  (total_differentiation_score + 0.01))
  )

# Developmental entropy
entropy_df <- prop_df %>%
  group_by(
    sample, response_code, response, response_group, response_binary,
    treatment, histology, stage, recurrence
  ) %>%
  summarise(
    developmental_entropy = -sum(prop * log(prop + 1e-9)),
    .groups = "drop"
  )

metric_df <- metric_df %>%
  left_join(
    entropy_df,
    by = c(
      "sample",
      "response_code",
      "response",
      "response_group",
      "response_binary",
      "treatment",
      "histology",
      "stage",
      "recurrence"
    )
  )

# Cycling blastema fraction
blastema_prolif_df <- obj.tr@meta.data %>%
  filter(developmental_group == "Blastema") %>%
  count(sample, Phase, name = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  filter(Phase %in% c("S", "G2M")) %>%
  group_by(sample) %>%
  summarise(
    cycling_blastema_fraction = sum(prop),
    .groups = "drop"
  )

metric_df <- metric_df %>%
  left_join(blastema_prolif_df, by = "sample") %>%
  mutate(
    cycling_blastema_fraction = ifelse(
      is.na(cycling_blastema_fraction),
      0,
      cycling_blastema_fraction
    )
  )

write.csv(
  metric_df,
  file.path(outdir, "plasticity_response_metrics_response123.csv"),
  row.names = FALSE
)

# ============================================================
# Statistics
# ============================================================

test_vars <- c(
  "blastema_persistence_score",
  "total_differentiation_score",
  "epithelial_differentiation_score",
  "stromal_myogenic_differentiation_score",
  "plasticity_log_ratio",
  "resistance_log_ratio",
  "developmental_entropy",
  "cycling_blastema_fraction"
)

stats_df <- lapply(test_vars, function(v) {
  
  x <- metric_df %>%
    filter(
      !is.na(.data[[v]]),
      response_group %in% c("Responder", "Non-responder")
    )
  
  if (length(unique(x$response_group)) < 2) {
    return(tibble(metric = v, p_value = NA_real_))
  }
  
  test <- wilcox.test(
    as.formula(paste(v, "~ response_group")),
    data = x,
    exact = FALSE
  )
  
  tibble(
    metric = v,
    p_value = test$p.value
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      is.na(p_adj) ~ NA_character_,
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001  ~ "***",
      p_adj < 0.01   ~ "**",
      p_adj < 0.05   ~ "*",
      TRUE ~ "ns"
    )
  )

write.csv(
  stats_df,
  file.path(outdir, "responder_vs_nonresponder_statistics_response123.csv"),
  row.names = FALSE
)

print(stats_df)

# ============================================================
# Colors
# ============================================================

resp_cols <- c(
  "Responder" = "#4DAF4A",
  "Non-responder" = "#E41A1C"
)

dev_cols <- c(
  "Blastema" = "#E64B35",
  "Epithelial" = "#4DBBD5",
  "Stromal/Myogenic" = "#00A087",
  "Other" = "grey80"
)

# ============================================================
# Plot 1: sample-wise composition
# ============================================================

sample_order <- metric_df %>%
  mutate(response_group = factor(response_group, levels = c("Non-responder", "Responder"))) %>%
  arrange(response_group, desc(blastema_persistence_score)) %>%
  pull(sample)

prop_df$sample <- factor(prop_df$sample, levels = sample_order)
prop_df$response_group <- factor(prop_df$response_group, levels = c("Non-responder", "Responder"))

p_sample_stack <- ggplot(
  prop_df,
  aes(x = sample, y = prop * 100, fill = developmental_group)
) +
  geom_col(color = "black", linewidth = 0.2) +
  facet_grid(. ~ response_group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = dev_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Malignant composition (%)",
    fill = "Developmental group",
    title = "Developmental-state composition in treated WT by tumor-burden response"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "01_samplewise_developmental_composition_response123.pdf"),
  p_sample_stack,
  width = 10,
  height = 5
)

# ============================================================
# Plot helper
# ============================================================

plot_metric <- function(df, metric, ylab, title) {
  
  pval <- stats_df %>%
    filter(.data$metric == metric) %>%
    pull(significance)
  
  ggplot(
    df,
    aes(x = response_group, y = .data[[metric]], fill = response_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.12, size = 3, alpha = 0.9) +
    scale_fill_manual(values = resp_cols) +
    labs(
      x = NULL,
      y = ylab,
      title = paste0(title, " (", pval, ")")
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11)
    )
}

p1 <- plot_metric(
  metric_df,
  "blastema_persistence_score",
  "Blastema fraction",
  "Blastema persistence"
)

p2 <- plot_metric(
  metric_df,
  "total_differentiation_score",
  "Differentiated fraction",
  "Total differentiation"
)

p3 <- plot_metric(
  metric_df,
  "plasticity_log_ratio",
  "log2((Differentiated + 0.01)/(Blastema + 0.01))",
  "Plasticity log-ratio"
)

p4 <- plot_metric(
  metric_df,
  "resistance_log_ratio",
  "log2((Blastema + 0.01)/(Differentiated + 0.01))",
  "Resistance log-ratio"
)

p5 <- plot_metric(
  metric_df,
  "developmental_entropy",
  "Developmental entropy",
  "State diversity"
)

p6 <- plot_metric(
  metric_df,
  "cycling_blastema_fraction",
  "Cycling blastema fraction",
  "Cycling blastema"
)

combined_metrics <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave(
  file.path(outdir, "02_response_plasticity_metrics_response123.pdf"),
  combined_metrics,
  width = 13,
  height = 8
)

# ============================================================
# Plot 3: fine state heatmap
# ============================================================

state_order <- c(
  "blastema 1",
  "S blastema 1",
  "mitotic blastema 1",
  "blastema 2",
  "mitotic blastema 2",
  "blastema 3",
  "tubules",
  "S tubules",
  "mitotic tubules",
  "podocyte-like",
  "PAX3+ myogenic precursors",
  "myoblast-like",
  "striated myocyte-like",
  "smooth myocyte-like precursors",
  "mitotic smooth myocyte-like precursors",
  "nascent smooth myocyte-like",
  "smooth myocyte-like",
  "mitotic smooth myocyte-like"
)

state_prop_df <- state_prop_df %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    response_group = factor(response_group, levels = c("Non-responder", "Responder")),
    state = factor(state, levels = rev(state_order))
  )

p_state_heat <- ggplot(
  state_prop_df,
  aes(x = sample, y = state, fill = prop)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(. ~ response_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "white",
    high = "black",
    labels = percent_format()
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Fraction",
    title = "Fine malignant cell-state composition by tumor-burden response"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "03_cell_state_composition_heatmap_response123.pdf"),
  p_state_heat,
  width = 12,
  height = 7
)

# ============================================================
# Plot 4: treatment-regimen stratified log-ratio
# ============================================================

metric_df$treatment <- factor(
  metric_df$treatment,
  levels = c("EE-4A", "DD-4A", "V+D", "I", "UNK")
)

p_treatment <- ggplot(
  metric_df,
  aes(x = treatment, y = plasticity_log_ratio, fill = response_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = response_group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    size = 3
  ) +
  scale_fill_manual(values = resp_cols, na.value = "grey80") +
  scale_color_manual(values = resp_cols, na.value = "grey80") +
  labs(
    x = "Treatment regimen",
    y = "Plasticity log-ratio",
    fill = "Response",
    color = "Response",
    title = "Developmental plasticity by regimen and tumor-burden response"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(outdir, "04_plasticity_log_ratio_by_treatment_response123.pdf"),
  p_treatment,
  width = 8,
  height = 5
)

# ============================================================
# Plot 5: resistance by histology
# ============================================================

p_histology <- ggplot(
  metric_df,
  aes(x = histology, y = resistance_log_ratio, fill = response_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = response_group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    size = 3
  ) +
  scale_fill_manual(values = resp_cols, na.value = "grey80") +
  scale_color_manual(values = resp_cols, na.value = "grey80") +
  labs(
    x = "Histology",
    y = "Resistance log-ratio",
    fill = "Response",
    color = "Response",
    title = "Blastema persistence and resistance by histology"
  ) +
  theme_classic(base_size = 13)

ggsave(
  file.path(outdir, "05_resistance_log_ratio_by_histology_response123.pdf"),
  p_histology,
  width = 7,
  height = 5
)

# ============================================================
# Plot 6: individual response categories 1/2/3
# ============================================================

metric_df$response <- factor(
  metric_df$response,
  levels = c(
    "Decreased tumor burden",
    "Unchanged tumor burden",
    "Increased tumor burden"
  )
)

p_response123 <- ggplot(
  metric_df,
  aes(x = response, y = blastema_persistence_score, fill = response)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 3) +
  labs(
    x = NULL,
    y = "Blastema fraction",
    title = "Blastema persistence across response categories"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "06_blastema_fraction_by_response_category_123.pdf"),
  p_response123,
  width = 7,
  height = 5
)

# ============================================================
# Optional UMAPs
# ============================================================

umap_reduction <- if ("umap" %in% names(obj@reductions)) {
  "umap"
} else if ("umap.recomputed" %in% names(obj@reductions)) {
  "umap.recomputed"
} else {
  NA_character_
}

if (!is.na(umap_reduction)) {
  
  p_umap_resp <- DimPlot(
    obj.tr,
    reduction = umap_reduction,
    group.by = "response_group",
    cols = resp_cols,
    pt.size = 0.15,
    raster = FALSE
  ) +
    theme_classic(base_size = 13) +
    ggtitle("Malignant cells by tumor-burden response")
  
  ggsave(
    file.path(outdir, "07_UMAP_response_group_response123.pdf"),
    p_umap_resp,
    width = 6,
    height = 5
  )
  
  p_umap_dev <- DimPlot(
    obj.tr,
    reduction = umap_reduction,
    group.by = "developmental_group",
    cols = dev_cols,
    pt.size = 0.15,
    raster = FALSE
  ) +
    theme_classic(base_size = 13) +
    ggtitle("Malignant developmental groups")
  
  ggsave(
    file.path(outdir, "08_UMAP_developmental_group_response123.pdf"),
    p_umap_dev,
    width = 6,
    height = 5
  )
}

# ============================================================
# Additional bar plots to test constrained plasticity
# ============================================================

# -----------------------------
# 1. Mean developmental composition by response group
# -----------------------------

group_bar_df <- prop_df %>%
  group_by(response_group, developmental_group) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sem = sd(prop, na.rm = TRUE) / sqrt(n_distinct(sample)),
    .groups = "drop"
  )

p_group_bar <- ggplot(
  group_bar_df,
  aes(x = response_group, y = mean_prop * 100, fill = developmental_group)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.75) +
  scale_fill_manual(values = dev_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Mean malignant composition (%)",
    fill = "Developmental group",
    title = "Mean developmental-state occupancy by tumor-burden response"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "09_group_mean_developmental_composition_response_bar.pdf"),
  p_group_bar,
  width = 6,
  height = 5
)

# -----------------------------
# 2. Early / intermediate / late developmental state grouping
# -----------------------------

early_states <- c(
  "blastema 1",
  "blastema 2",
  "blastema 3",
  "S blastema 1",
  "mitotic blastema 1",
  "mitotic blastema 2"
)

intermediate_states <- c(
  "tubules",
  "mitotic tubules",
  "podocyte-like",
  "PAX3+ myogenic precursors",
  "myoblast-like"
)

late_states <- c(
  "smooth myocyte-like",
  "smooth myocyte-like precursors",
  "nascent smooth myocyte-like",
  "mitotic smooth myocyte-like",
  "mitotic smooth myocyte-like precursors",
  "striated myocyte-like"
)

obj.tr$developmental_stage <- case_when(
  obj.tr$state %in% early_states ~ "Early blastemal",
  obj.tr$state %in% intermediate_states ~ "Intermediate differentiation",
  obj.tr$state %in% late_states ~ "Late stromal/myogenic",
  TRUE ~ "Other"
)

obj.tr$developmental_stage <- factor(
  obj.tr$developmental_stage,
  levels = c(
    "Early blastemal",
    "Intermediate differentiation",
    "Late stromal/myogenic",
    "Other"
  )
)

stage_prop_df <- obj.tr@meta.data %>%
  count(
    sample,
    response_group,
    response,
    treatment,
    histology,
    developmental_stage,
    name = "n"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_malignant = sum(n),
    prop = n / total_malignant
  ) %>%
  ungroup()

stage_cols <- c(
  "Early blastemal" = "#D73027",
  "Intermediate differentiation" = "#FDAE61",
  "Late stromal/myogenic" = "#1A9850",
  "Other" = "grey80"
)

stage_bar_df <- stage_prop_df %>%
  group_by(response_group, developmental_stage) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  )

p_stage_bar <- ggplot(
  stage_bar_df,
  aes(x = response_group, y = mean_prop * 100, fill = developmental_stage)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.75) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Mean malignant composition (%)",
    fill = "Developmental stage",
    title = "Early versus late developmental-state occupancy"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "10_early_intermediate_late_response_bar.pdf"),
  p_stage_bar,
  width = 7,
  height = 5
)

# -----------------------------
# 3. Same early / late plot by histology
# -----------------------------

stage_bar_histology_df <- stage_prop_df %>%
  group_by(histology, developmental_stage) %>%
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  )

p_stage_histology <- ggplot(
  stage_bar_histology_df,
  aes(x = histology, y = mean_prop * 100, fill = developmental_stage)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.75) +
  scale_fill_manual(values = stage_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Mean malignant composition (%)",
    fill = "Developmental stage",
    title = "Developmental constraint by histology"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "11_early_intermediate_late_histology_bar.pdf"),
  p_stage_histology,
  width = 6,
  height = 5
)

# -----------------------------
# 4. Early vs late progression score
# -----------------------------

progression_df <- stage_prop_df %>%
  filter(developmental_stage %in% c(
    "Early blastemal",
    "Intermediate differentiation",
    "Late stromal/myogenic"
  )) %>%
  select(sample, response_group, response, treatment, histology,
         developmental_stage, prop) %>%
  pivot_wider(
    names_from = developmental_stage,
    values_from = prop,
    values_fill = 0
  ) %>%
  mutate(
    developmental_progression_score =
      `Late stromal/myogenic` + `Intermediate differentiation` - `Early blastemal`,
    constrained_plasticity_score =
      `Early blastemal` - (`Intermediate differentiation` + `Late stromal/myogenic`)
  )

write.csv(
  progression_df,
  file.path(outdir, "developmental_progression_constraint_scores.csv"),
  row.names = FALSE
)

p_progression <- ggplot(
  progression_df,
  aes(x = response_group, y = developmental_progression_score, fill = response_group)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.12, size = 3) +
  scale_fill_manual(values = resp_cols) +
  labs(
    x = NULL,
    y = "Progression score\n(intermediate + late − early)",
    title = "Developmental progression by response"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "12_developmental_progression_score_response.pdf"),
  p_progression,
  width = 5,
  height = 5
)

p_constraint_histology <- ggplot(
  progression_df,
  aes(x = histology, y = constrained_plasticity_score, fill = histology)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.12, size = 3, aes(shape = response_group)) +
  labs(
    x = NULL,
    y = "Constraint score\n(early − intermediate − late)",
    title = "Constrained plasticity by histology"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(outdir, "13_constrained_plasticity_score_histology.pdf"),
  p_constraint_histology,
  width = 5.5,
  height = 5
)

message("Additional constrained-plasticity bar plots generated.")

# ============================================================
# Save final object
# ============================================================

saveRDS(
  obj,
  file.path(outdir, "Wilms_TME_with_response123_plasticity_metadata.rds")
)

message("Done.")
message("Outputs written to: ", outdir)
message("Main files:")
message("  01_samplewise_developmental_composition_response123.pdf")
message("  02_response_plasticity_metrics_response123.pdf")
message("  03_cell_state_composition_heatmap_response123.pdf")
message("  04_plasticity_log_ratio_by_treatment_response123.pdf")
message("  05_resistance_log_ratio_by_histology_response123.pdf")
message("  plasticity_response_metrics_response123.csv")
