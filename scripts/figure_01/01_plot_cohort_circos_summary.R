#!/usr/bin/env Rscript

# =============================================================================
# Script: 01_plot_cohort_circos_summary.R
# Project: Wilms tumor plasticity manuscript
# Figure: Figure 1 - cohort-level circos summary
#
# Purpose:
#   Generate a manuscript-ready circos plot summarizing patient-level clinical,
#   pathology, genomic, and malignant-compartment composition features.
#
# Inputs:
#   data/metadata/cohort.tsv
#       Patient-level metadata. Required columns include:
#       patient, race, sex, age, treatment, tumor, survival, stage, recurrence,
#       laterality, size, weight, histology, focality, neoadjuvant, response.
#
#   data/metadata/FGA_TMB.tsv
#       Patient-level genomic summary. Required columns:
#       Patient, FGA, TMB.
#
#   data/processed/TME.RData
#       Seurat object named `tme_`. Required metadata columns:
#       orig.ident, subcompartment.
#
# Outputs:
#   figures/figure_01/cohort_summary_circos.pdf
#   figures/figure_01/cohort_summary_circos_legend_only.pdf
#
# Notes:
#   - Paths are relative to the repository root.
#   - Run this script from the repository root:
#       Rscript scripts/figure_01/01_plot_cohort_circos_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(gridBase)
  library(tidyverse)
  library(Seurat)
})

# -----------------------------------------------------------------------------
# User-configurable paths
# -----------------------------------------------------------------------------
cohort_path <- file.path("data", "metadata", "cohort.tsv")
fga_tmb_path <- file.path("data", "metadata", "FGA_TMB.tsv")
tme_path <- file.path("data", "processed", "TME.RData")
out_dir <- file.path("figures", "figure_01")

out_pdf <- file.path(out_dir, "cohort_summary_circos.pdf")
legend_pdf <- file.path(out_dir, "cohort_summary_circos_legend_only.pdf")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input file: ", path, call. = FALSE)
  }
}

safe_cat_cols <- function(values, cols, na_col = "white") {
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "unk"

  cols <- c(cols, "unk" = na_col)
  bg_cols <- cols[values]
  bg_cols[is.na(bg_cols)] <- na_col
  bg_cols
}

add_cat_track <- function(values, cols, sectors, border = "black", height = 0.014) {
  bg_cols <- safe_cat_cols(values, cols)
  circos.trackPlotRegion(
    sectors = sectors,
    ylim = c(0, 1),
    track.height = height,
    bg.col = bg_cols,
    bg.border = border,
    panel.fun = function(x, y) NULL
  )
}

add_cont_track <- function(values, col_fun, sectors, border = "black", height = 0.014) {
  bg_cols <- col_fun(values)
  bg_cols[is.na(values)] <- "white"

  circos.trackPlotRegion(
    sectors = sectors,
    ylim = range(values, na.rm = TRUE),
    track.height = height,
    bg.col = bg_cols,
    bg.border = border,
    panel.fun = function(x, y) NULL
  )
}

make_cat_legend <- function(labels, cols, title, nrow = 1, ncol = NULL,
                            border = NA, background = NULL) {
  Legend(
    at = labels,
    legend_gp = gpar(col = cols),
    type = "points",
    pch = 15,
    size = unit(3, "mm"),
    labels_gp = gpar(fontsize = 5.5),
    title_gp = gpar(fontsize = 6.5, fontface = "bold"),
    title = title,
    title_position = "topleft",
    background = background,
    border = border,
    nrow = nrow,
    ncol = ncol
  )
}

make_cont_legend <- function(values, col_fun, title) {
  Legend(
    at = round(c(min(values, na.rm = TRUE), mean(range(values, na.rm = TRUE)), max(values, na.rm = TRUE)), 2),
    col_fun = col_fun,
    title_position = "topleft",
    title = title,
    direction = "horizontal",
    legend_width = unit(18, "mm"),
    labels_gp = gpar(fontsize = 5.5),
    title_gp = gpar(fontsize = 6.5, fontface = "bold")
  )
}

# -----------------------------------------------------------------------------
# Validate and load inputs
# -----------------------------------------------------------------------------
check_file_exists(cohort_path)
check_file_exists(fga_tmb_path)
check_file_exists(tme_path)

cohort <- read.delim(
  file = cohort_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  as.is = TRUE,
  check.names = FALSE
)

fg_tmb <- read.delim(
  file = fga_tmb_path,
  header = TRUE,
  sep = "\t",
  quote = "",
  as.is = TRUE,
  check.names = FALSE
)

required_cohort_cols <- c(
  "patient", "race", "sex", "age", "treatment", "tumor", "survival",
  "stage", "recurrence", "laterality", "size", "weight", "histology",
  "focality", "neoadjuvant", "response"
)
required_fga_cols <- c("Patient", "FGA", "TMB")
missing_cohort_cols <- setdiff(required_cohort_cols, colnames(cohort))
missing_fga_cols <- setdiff(required_fga_cols, colnames(fg_tmb))

if (length(missing_cohort_cols) > 0) {
  stop("Missing columns in cohort.tsv: ", paste(missing_cohort_cols, collapse = ", "), call. = FALSE)
}
if (length(missing_fga_cols) > 0) {
  stop("Missing columns in FGA_TMB.tsv: ", paste(missing_fga_cols, collapse = ", "), call. = FALSE)
}

rownames(cohort) <- cohort$patient
rownames(fg_tmb) <- fg_tmb$Patient

cohort$FGA <- fg_tmb[cohort$patient, "FGA"]
cohort$TMB <- fg_tmb[cohort$patient, "TMB"]

cohort$FGA <- as.numeric(cohort$FGA)
cohort$TMB <- as.numeric(cohort$TMB)
cohort$age <- as.numeric(cohort$age)
cohort$stage <- as.numeric(cohort$stage)
cohort$weight <- as.numeric(cohort$weight)

# Untreated samples are not assigned a neoadjuvant regimen.
cohort$neoadjuvant[cohort$treatment %in% "untreated"] <- "unk"

load(tme_path)
if (!exists("tme_")) {
  stop("The TME.RData file must contain a Seurat object named `tme_`.", call. = FALSE)
}

tme_meta <- tme_@meta.data
required_tme_cols <- c("orig.ident", "subcompartment")
missing_tme_cols <- setdiff(required_tme_cols, colnames(tme_meta))
if (length(missing_tme_cols) > 0) {
  stop("Missing columns in tme_@meta.data: ", paste(missing_tme_cols, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# Compute malignant subcompartment percentages per sample
# -----------------------------------------------------------------------------
tme_df <- tme_meta %>%
  dplyr::filter(subcompartment %in% c("blastema", "stroma", "epithelium")) %>%
  dplyr::count(orig.ident, subcompartment, name = "n") %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(frequency = round(100 * n / sum(n), 2)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(patient = orig.ident) %>%
  tidyr::complete(
    patient = cohort$patient,
    subcompartment = c("blastema", "stroma", "epithelium"),
    fill = list(n = 0, frequency = 0)
  )

tme_df$subcompartment <- factor(
  tme_df$subcompartment,
  levels = c("blastema", "stroma", "epithelium")
)

# -----------------------------------------------------------------------------
# Color palettes
# -----------------------------------------------------------------------------
race_cols <- c(
  "Black,Hispanic/Latino" = "yellow4",
  "White" = "pink",
  "Black" = "violet",
  "Asian" = "purple2",
  "Hispanic/Latino" = "purple4"
)

sex_cols <- c("female" = "lightgreen", "male" = "darkgreen")
treatment_cols <- c("treated" = "red", "untreated" = "blue")
tumor_cols <- c("Primary" = "#6BAED6", "Recurrence" = "#FB6A4A")
survival_cols <- c("dead" = "black", "alive" = "grey75")
recurrence_cols <- c("yes" = "black", "no" = "yellow")
laterality_cols <- c("left" = "#D2B48C", "right" = "skyblue", "bilateral" = "#4A2F26")
size_cols <- c("<=10" = "orange", ">=10" = "black")
histology_cols <- c("favorable" = "green", "anaplastic" = "black")
focality_cols <- c("unifocal" = "#602CA1", "multifocal" = "#BF5842")
neoadj_cols <- c("EE-4A" = "orange", "DD-4A" = "red", "Vincristine+Doxorubicin" = "#AA9922", "I" = "grey35")
response_cols <- c("decreased" = "green2", "increased" = "red", "unchanged" = "grey50")

col_fun_tme <- list(
  blastema = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == "blastema"], na.rm = TRUE), c("pink", "red3")),
  stroma = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == "stroma"], na.rm = TRUE), c("lightgreen", "green4")),
  epithelium = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == "epithelium"], na.rm = TRUE), c("lightblue", "blue3"))
)

col_fun_fga <- colorRamp2(range(cohort$FGA, na.rm = TRUE), c("pink", "purple4"))
col_fun_tmb <- colorRamp2(range(cohort$TMB, na.rm = TRUE), c("khaki", "brown4"))
col_fun_age <- colorRamp2(range(cohort$age, na.rm = TRUE), c("yellow", "darkorange4"))
col_fun_stage <- colorRamp2(range(cohort$stage, na.rm = TRUE), c("pink", "purple"))
col_fun_weight <- colorRamp2(range(cohort$weight, na.rm = TRUE), c("lightblue", "darkblue"))

# -----------------------------------------------------------------------------
# Draw circos plot
# -----------------------------------------------------------------------------
graphics.off()
circos.clear()

pdf(file = out_pdf, width = 15.71, height = 14.81)
par(mar = c(1, 1, 1, 1))

sector_ids <- cohort$patient

circos.par(
  start.degree = 90,
  gap.degree = c(rep(1, length(sector_ids) - 1), 90),
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.001, 0.001),
  canvas.xlim = c(-1.25, 1.25),
  canvas.ylim = c(-1.25, 1.25)
)

circos.initialize(
  factors = sector_ids,
  xlim = cbind(rep(0, length(sector_ids)), rep(3, length(sector_ids)))
)

# Tumor subcompartment cellularity bar track.
circos.trackPlotRegion(
  sectors = sector_ids,
  ylim = c(0, max(tme_df$frequency, na.rm = TRUE) * 1.2),
  track.height = 0.12,
  bg.border = "black",
  bg.col = NA,
  panel.fun = function(x, y) {
    sector_index <- CELL_META$sector.index
    patient_data <- tme_df[tme_df$patient == sector_index, ]
    bar_positions <- c(0.5, 1.5, 2.5)

    for (i in seq_along(bar_positions)) {
      subcompartment_i <- levels(tme_df$subcompartment)[i]
      feature_data <- patient_data[patient_data$subcompartment == subcompartment_i, ]

      if (nrow(feature_data) == 1 && !is.na(feature_data$frequency)) {
        circos.barplot(
          value = feature_data$frequency,
          bar_width = 0.85,
          pos = bar_positions[i],
          col = col_fun_tme[[subcompartment_i]](feature_data$frequency),
          border = NA
        )
      }
    }
  }
)

# Genomic burden tracks.
for (feature in c("FGA", "TMB")) {
  col_fun <- if (feature == "FGA") col_fun_fga else col_fun_tmb
  circos.trackPlotRegion(
    sectors = sector_ids,
    ylim = c(0, max(cohort[[feature]], na.rm = TRUE) * 1.2),
    track.height = 0.055,
    bg.border = "black",
    bg.col = NA,
    panel.fun = function(x, y) {
      val <- cohort[[feature]][CELL_META$sector.numeric.index]
      if (!is.na(val)) {
        circos.barplot(value = val, bar_width = 0.85, pos = 1.5, col = col_fun(val), border = NA)
      }
    }
  )
}

# Annotation tracks.
add_cat_track(cohort$race, race_cols, sector_ids, border = NA)
add_cat_track(cohort$sex, sex_cols, sector_ids, border = NA)
add_cont_track(cohort$age, col_fun_age, sector_ids, border = NA)

treatment_bg <- safe_cat_cols(cohort$treatment, treatment_cols)
font_cols <- ifelse(cohort$treatment == "treated", "white", "black")

circos.trackPlotRegion(
  sectors = sector_ids,
  ylim = c(0, 1),
  track.height = 0.045,
  bg.col = treatment_bg,
  bg.border = NA,
  panel.fun = function(x, y) {
    i <- CELL_META$sector.numeric.index
    circos.text(
      x = CELL_META$xcenter,
      y = 0.5,
      labels = cohort$patient[i],
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0.5, 0.5),
      cex = 0.75,
      col = font_cols[i],
      font = 2
    )
  }
)

add_cat_track(cohort$tumor, tumor_cols, sector_ids, border = NA, height = 0.018)
add_cat_track(cohort$survival, survival_cols, sector_ids)
add_cont_track(cohort$stage, col_fun_stage, sector_ids)
add_cat_track(cohort$recurrence, recurrence_cols, sector_ids)
add_cat_track(cohort$laterality, laterality_cols, sector_ids)
add_cat_track(cohort$size, size_cols, sector_ids)
add_cont_track(cohort$weight, col_fun_weight, sector_ids)
add_cat_track(cohort$histology, histology_cols, sector_ids)
add_cat_track(cohort$focality, focality_cols, sector_ids, border = "white")
add_cat_track(cohort$neoadjuvant, neoadj_cols, sector_ids)
add_cat_track(cohort$response, response_cols, sector_ids)

# -----------------------------------------------------------------------------
# Legends
# -----------------------------------------------------------------------------
blastema_lgd <- make_cont_legend(tme_df$frequency[tme_df$subcompartment == "blastema"], col_fun_tme$blastema, "Blastema content (%)")
stroma_lgd <- make_cont_legend(tme_df$frequency[tme_df$subcompartment == "stroma"], col_fun_tme$stroma, "Stroma content (%)")
epithelium_lgd <- make_cont_legend(tme_df$frequency[tme_df$subcompartment == "epithelium"], col_fun_tme$epithelium, "Epithelium content (%)")
fga_lgd <- make_cont_legend(cohort$FGA, col_fun_fga, "FGA (%)")
tmb_lgd <- make_cont_legend(cohort$TMB, col_fun_tmb, "TMB (Mut/Mb)")
age_lgd <- make_cont_legend(cohort$age, col_fun_age, "Age (months)")
stage_lgd <- make_cont_legend(cohort$stage, col_fun_stage, "Cancer stage")
weight_lgd <- make_cont_legend(cohort$weight, col_fun_weight, "Tumor weight (g)")

race_lgd <- make_cat_legend(names(race_cols), unname(race_cols), "Race", nrow = NULL, ncol = 1)
sex_lgd <- make_cat_legend(c("Female", "Male"), unname(sex_cols[c("female", "male")]), "Sex")
treatment_lgd <- make_cat_legend(c("Treated", "Untreated"), unname(treatment_cols[c("treated", "untreated")]), "Patient treatment", border = "black")
tumor_lgd <- make_cat_legend(c("Primary", "Recurrence"), unname(tumor_cols[c("Primary", "Recurrence")]), "Tumor", border = "black")
survival_lgd <- make_cat_legend(c("Dead", "Alive"), unname(survival_cols[c("dead", "alive")]), "Survival")
recurrence_lgd <- make_cat_legend(c("Yes", "No"), unname(recurrence_cols[c("yes", "no")]), "Recurrence")
laterality_lgd <- make_cat_legend(c("Left", "Right", "Bilateral"), unname(laterality_cols[c("left", "right", "bilateral")]), "Laterality")
size_lgd <- make_cat_legend(c("<=10", ">=10"), unname(size_cols[c("<=10", ">=10")]), "Tumor size (cm)")
histology_lgd <- make_cat_legend(c("Favorable", "Anaplastic"), unname(histology_cols[c("favorable", "anaplastic")]), "Histology")
focality_lgd <- make_cat_legend(c("Unifocal", "Multifocal"), unname(focality_cols[c("unifocal", "multifocal")]), "Tumor focality")
neoadj_lgd <- make_cat_legend(names(neoadj_cols), unname(neoadj_cols), "Neoadjuvant", nrow = NULL, ncol = 1, border = "black")
response_lgd <- make_cat_legend(names(response_cols), unname(response_cols), "Response", nrow = NULL, ncol = 1, border = "black")

lgd_list_left <- packLegend(
  direction = "vertical", gap = unit(1.5, "mm"),
  race_lgd, sex_lgd, age_lgd, treatment_lgd, tumor_lgd,
  survival_lgd, recurrence_lgd, laterality_lgd, size_lgd,
  histology_lgd, focality_lgd
)

lgd_list_right <- packLegend(
  direction = "vertical", gap = unit(1.5, "mm"),
  stage_lgd, weight_lgd, neoadj_lgd, response_lgd,
  blastema_lgd, stroma_lgd, epithelium_lgd, fga_lgd, tmb_lgd
)

# Overlay grid legends onto the base/circos graphics device.
vps <- gridBase::baseViewports()
pushViewport(vps$inner)
pushViewport(vps$figure)
pushViewport(vps$plot)

draw(lgd_list_left, x = unit(0.02, "npc"), y = unit(0.83, "npc"), just = c("left", "top"))
draw(lgd_list_right, x = unit(0.22, "npc"), y = unit(0.83, "npc"), just = c("left", "top"))

upViewport(3)

circos.clear()
dev.off()
graphics.off()

# -----------------------------------------------------------------------------
# Separate legend-only PDF
# -----------------------------------------------------------------------------
pdf(legend_pdf, width = 8.5, height = 11)
grid.newpage()

lgd_list_all <- packLegend(
  direction = "vertical", gap = unit(2, "mm"),
  blastema_lgd, stroma_lgd, epithelium_lgd, fga_lgd, tmb_lgd,
  race_lgd, sex_lgd, age_lgd, treatment_lgd, tumor_lgd,
  survival_lgd, stage_lgd, recurrence_lgd, laterality_lgd, size_lgd,
  weight_lgd, histology_lgd, focality_lgd, neoadj_lgd, response_lgd
)

draw(lgd_list_all, x = unit(0.05, "npc"), y = unit(0.98, "npc"), just = c("left", "top"))
dev.off()

message("Saved circos plot: ", out_pdf)
message("Saved legend-only PDF: ", legend_pdf)
