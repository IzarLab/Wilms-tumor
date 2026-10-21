# This script plots a circos summary of the cohort: patient-level clinical, pathology, genomic (FGA, TMB) and tumor subcompartment
# (blastema, stroma, epithelium as % of all TME cells) features; each track is a feature and each sector a patient
# Run from the root of the repository: Rscript scripts/Fig_1/b.R

library(circlize)
library(ComplexHeatmap)
library(grid)
library(tidyverse)
library(Seurat)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in data ####

## step 1.1: cohort metadata, FGA and TMB ####

cohort = read.delim(file = 'data/Misc/cohort.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
fg_tmb = read.delim(file = 'data/Misc/FGA_TMB.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
rownames(cohort) = cohort$patient
rownames(fg_tmb) = fg_tmb$Patient
cohort$FGA = as.numeric(fg_tmb[cohort$patient, 'FGA'])     # patients absent from FGA_TMB.tsv get NA
cohort$TMB = as.numeric(fg_tmb[cohort$patient, 'TMB'])
cohort$weight = suppressWarnings(as.numeric(cohort$weight))     # 'unk' weights become NA
rm(fg_tmb)

## step 1.2: TME metadata ####

load('data/compartments/TME.RData')     # Seurat object tme_
tme_meta = tme_@meta.data
rm(tme_); gc()

# STEP 2: Tumor subcompartment percentages per sample (of all TME cells) ####

subcompartments = c('blastema', 'stroma', 'epithelium')
tme_df = tme_meta %>%
         dplyr::filter(orig.ident %in% cohort$patient) %>%
         dplyr::group_by(patient = orig.ident) %>%
         dplyr::summarise(blastema = round(100*sum(subcompartment %in% 'blastema')/dplyr::n(), 2),
                          stroma = round(100*sum(subcompartment %in% 'stroma')/dplyr::n(), 2),
                          epithelium = round(100*sum(subcompartment %in% 'epithelium')/dplyr::n(), 2),
                          .groups = 'drop') %>%
         tidyr::pivot_longer(cols = dplyr::all_of(subcompartments), names_to = 'subcompartment', values_to = 'frequency') %>%
         as.data.frame()
tme_df$subcompartment = factor(tme_df$subcompartment, levels = subcompartments)
rm(tme_meta); gc()

# STEP 3: Plotting the circos ####

## step 3.1: colors ####

race_cols = c('Black,Hispanic/Latino' = 'yellow4', 'White' = 'pink', 'Black' = 'violet', 'Asian' = 'purple2', 'Hispanic/Latino' = 'purple4')
sex_cols = c('female' = 'lightgreen', 'male' = 'darkgreen')
treatment_cols = c('treated' = 'red', 'untreated' = 'blue')
survival_cols = c('dead' = 'black', 'alive' = 'grey75')
recurrence_cols = c('yes' = 'black', 'no' = 'yellow')
tumor_cols = c('1' = 'blue', '2' = 'orange')
laterality_cols = c('left' = '#D2B48C', 'right' = 'skyblue', 'bilateral' = '#4A2F26')
size_cols = c('<=10' = 'orange', '>=10' = 'black')
histology_cols = c('favorable' = 'green', 'anaplastic' = 'black')
focality_cols = c('unifocal' = '#602CA1', 'multifocal' = '#BF5842')
neoadj_cols = c('EE-4A' = 'orange', 'DD-4A' = 'red', 'Vincristine+Doxorubicin' = '#AA9922', 'I' = 'grey35')
response_cols = c('decreased' = 'green2', 'increased' = 'red', 'unchanged' = 'grey50')

lookup_cols = function(values, cols){ unname(cols[as.character(values)]) }     # values without a color (e.g. 'unk', NA) are left unfilled

col_fun_tme = list(blastema = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == 'blastema'], na.rm = T), c('pink', 'red3')),
                   stroma = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == 'stroma'], na.rm = T), c('lightgreen', 'green4')),
                   epithelium = colorRamp2(range(tme_df$frequency[tme_df$subcompartment == 'epithelium'], na.rm = T), c('lightblue', 'blue3')))
col_fun_fga = colorRamp2(range(cohort$FGA, na.rm = T), c('pink', 'purple4'))
col_fun_tmb = colorRamp2(range(cohort$TMB, na.rm = T), c('khaki', 'brown4'))
col_fun_age = colorRamp2(range(cohort$age, na.rm = T), c('yellow', 'darkorange4'))
col_fun_stage = colorRamp2(range(cohort$stage, na.rm = T), c('pink', 'purple'))
col_fun_weight = colorRamp2(range(cohort$weight, na.rm = T), c('lightblue', 'darkblue'))

## step 3.2: track helpers (used by many tracks) ####

add_annotation_track = function(bg_cols, sectors, border = 'black', height = 0.026)     # track filled with one background color per sector
{
  circos.trackPlotRegion(sectors = sectors, ylim = c(0, 1), track.height = height, bg.col = bg_cols, bg.border = border, panel.fun = function(x, y){ NULL })
}

add_bar_axis = function(values, n_sectors)     # dotted grid lines and y-axis ticks for bar tracks; tick labels on the first sector only
{
  breaks_ = round(seq(min(values, na.rm = T), max(values, na.rm = T), length.out = 3), 2)
  i_ = CELL_META$sector.numeric.index
  for(b_ in breaks_)
  {
    circos.lines(x = CELL_META$cell.xlim, y = rep(b_, 2), lty = 3, col = 'black')
    if(i_ == 1){ circos.yaxis(side = 'left', at = b_, labels = T) }
    if(1 < i_ & i_ < n_sectors){ circos.yaxis(side = 'left', at = b_, labels = F) }
    if(i_ == n_sectors){ circos.yaxis(side = 'right', at = b_, labels = F) }
  }
}

## step 3.3: circos layout ####

graphics.off()
circos.clear()
pdf(file = file.path(outdir_, 'cohort_summary_circos.pdf'), width = 15.71, height = 14.81)
par(mar = c(1, 1, 1, 1))

sector_ids = cohort$patient
n_sectors = length(sector_ids)
circos.par(start.degree = 90,
           gap.degree = c(rep(1, n_sectors - 1), 90),     # 90-degree gap for the legends
           cell.padding = c(0, 0, 0, 0),
           track.height = 0.1,                            # default track height
           track.margin = c(0.0, 0.008))
circos.initialize(factors = sector_ids, xlim = c(0, 3))

## step 3.4: genomic tracks ####

# 1a: tumor subcompartment content (blastema, stroma, epithelium bars)
circos.trackPlotRegion(sectors = sector_ids, ylim = c(0, max(tme_df$frequency)*1.2), track.height = 0.15, bg.border = 'black', bg.col = NA,
                       panel.fun = function(x, y)
                       {
                         patient_data = tme_df[tme_df$patient == CELL_META$sector.index, ]
                         bar_positions = c(0.5, 1.5, 2.5)
                         for(i_ in seq_along(bar_positions))
                         {
                           val_ = patient_data$frequency[patient_data$subcompartment == subcompartments[i_]]
                           circos.barplot(value = val_, bar_width = 0.85, pos = bar_positions[i_], col = col_fun_tme[[subcompartments[i_]]](val_), border = NA)
                           add_bar_axis(tme_df$frequency, n_sectors)     # redrawn after each bar to match the original figure's line weight
                         }
                       })

# 1b: FGA, 1c: TMB
for(feature_ in c('FGA', 'TMB'))
{
  local({
    values_ = cohort[[feature_]]
    col_fun = if(feature_ == 'FGA'){ col_fun_fga }else{ col_fun_tmb }
    circos.trackPlotRegion(sectors = sector_ids, ylim = c(0, max(values_, na.rm = T)*1.2), bg.border = 'black', bg.col = NA,
                           panel.fun = function(x, y)
                           {
                             val_ = values_[CELL_META$sector.numeric.index]
                             circos.barplot(value = val_, bar_width = 0.85, pos = 1.5, col = col_fun(val_), border = NA)
                             add_bar_axis(values_, n_sectors)
                           })
  })
}

## step 3.5: demographic tracks ####

add_annotation_track(lookup_cols(cohort$race, race_cols), sector_ids, border = NA, height = 0.027)     # 2: race
add_annotation_track(lookup_cols(cohort$sex, sex_cols), sector_ids, border = NA, height = 0.027)       # 3: sex
add_annotation_track(col_fun_age(cohort$age), sector_ids, border = NA, height = 0.027)                 # 4: age

# 5: patient ID, colored by treatment
font_cols = ifelse(cohort$treatment == 'treated', 'black', 'white')
circos.trackPlotRegion(sectors = sector_ids, ylim = c(0, 1), bg.col = lookup_cols(cohort$treatment, treatment_cols), bg.border = NA,
                       panel.fun = function(x, y)
                       {
                         i_ = CELL_META$sector.numeric.index
                         circos.text(x = CELL_META$xcenter, y = 0.5, labels = cohort$patient[i_], facing = 'clockwise', niceFacing = T,
                                     adj = c(0.5, 0.5), cex = 1, col = font_cols[i_], font = 2)
                       })

## step 3.6: pathology tracks ####

add_annotation_track(lookup_cols(cohort$survival, survival_cols), sector_ids)                   # 6: survival
add_annotation_track(col_fun_stage(cohort$stage), sector_ids)                                   # 7: cancer stage (I-V)
add_annotation_track(lookup_cols(cohort$recurrence, recurrence_cols), sector_ids)               # 8: recurrence
add_annotation_track(lookup_cols(cohort$tumor, tumor_cols), sector_ids)                         # 9: tumor (1 = primary, 2 = recurrence)
add_annotation_track(lookup_cols(cohort$laterality, laterality_cols), sector_ids)               # 10: tumor laterality
add_annotation_track(lookup_cols(cohort$size, size_cols), sector_ids)                           # 11: tumor size (<=10, >=10)
add_annotation_track(col_fun_weight(cohort$weight), sector_ids)                                 # 12: tumor weight
add_annotation_track(lookup_cols(cohort$histology, histology_cols), sector_ids)                 # 13: tumor histology (favorable, anaplastic)
add_annotation_track(lookup_cols(cohort$focality, focality_cols), sector_ids, border = 'white') # 14: tumor focality (unifocal, multifocal)
add_annotation_track(lookup_cols(cohort$neoadjuvant, neoadj_cols), sector_ids)                  # 15: neoadjuvant
add_annotation_track(lookup_cols(cohort$response, response_cols), sector_ids)                   # 16: response to neoadjuvant

## step 3.7: legends ####

make_cat_legend = function(labels, cols, title, nrow = 1, ncol = 1, border = NULL, background = NULL)     # categorical legend
{
  Legend(at = labels, legend_gp = gpar(col = cols), type = 'points', pch = 15, size = unit(5, 'mm'), title = title, title_position = 'topleft',
         background = background, border = border, nrow = nrow, ncol = ncol)
}
make_cont_legend = function(at, col_fun, title, border = NULL)     # continuous legend
{
  Legend(at = at, col_fun = col_fun, title_position = 'topleft', title = title, border = border, direction = 'horizontal')
}
min_mid_max = function(values){ c(min(values, na.rm = T), sum(range(values, na.rm = T))/2, max(values, na.rm = T)) }

# genomic events
blastema_lgd = make_cont_legend(round(min_mid_max(tme_df$frequency[tme_df$subcompartment == 'blastema']), 2), col_fun_tme$blastema, 'Blastema content (%)')
stroma_lgd = make_cont_legend(round(min_mid_max(tme_df$frequency[tme_df$subcompartment == 'stroma']), 2), col_fun_tme$stroma, 'Stroma content (%)')
epithelium_lgd = make_cont_legend(round(min_mid_max(tme_df$frequency[tme_df$subcompartment == 'epithelium']), 2), col_fun_tme$epithelium, 'Epithelium content (%)')
fga_lgd = make_cont_legend(round(min_mid_max(cohort$FGA), 2), col_fun_fga, 'FGA (%)')
tmb_lgd = make_cont_legend(round(min_mid_max(cohort$TMB), 2), col_fun_tmb, 'TMB (Mut/Mb)')

# demographics
race_lgd = make_cat_legend(names(race_cols), unname(race_cols), 'Race', nrow = NULL, border = NA)
sex_lgd = make_cat_legend(c('Female', 'Male'), unname(sex_cols), 'Sex', border = NA)
age_lgd = make_cont_legend(round(min_mid_max(cohort$age)), col_fun_age, 'Age (months)', border = NA)
treatment_lgd = make_cat_legend(c('Treated', 'Untreated'), unname(treatment_cols), 'Patient', border = 'black')

# pathology
survival_lgd = make_cat_legend(c('dead', 'alive'), unname(survival_cols), 'Survival')
stage_lgd = make_cont_legend(unique(sort(cohort$stage)), col_fun_stage, 'Cancer stage')
recurrence_lgd = make_cat_legend(c('Yes', 'No'), unname(recurrence_cols), 'Recurrence', background = 'black')
tumor_lgd = make_cat_legend(c('Primary', 'Recurrence'), unname(tumor_cols), 'Tumor')
laterality_lgd = make_cat_legend(c('left', 'right', 'bilateral'), unname(laterality_cols), 'Laterality')
size_lgd = make_cat_legend(c('<=10', '>=10'), unname(size_cols), 'Tumor size (cm)')
weight_lgd = make_cont_legend(min_mid_max(cohort$weight), col_fun_weight, 'Tumor weight (g)')
histology_lgd = make_cat_legend(c('Favorable', 'Anaplastic'), unname(histology_cols), 'Histology')
focality_lgd = make_cat_legend(c('Unifocal', 'Multifocal'), unname(focality_cols), 'Tumor focality', background = 'black')
neoadj_lgd = make_cat_legend(names(neoadj_cols), unname(neoadj_cols), 'Neoadjuvant', nrow = NULL, border = 'black')
response_lgd = make_cat_legend(names(response_cols), unname(response_cols), 'Response', nrow = NULL, border = 'black')

draw(packLegend(direction = 'vertical', survival_lgd, stage_lgd, recurrence_lgd, tumor_lgd, laterality_lgd,
                size_lgd, weight_lgd, histology_lgd, focality_lgd, neoadj_lgd, response_lgd), x = unit(30, 'mm'), y = unit(360, 'mm'), just = c('top'))     # 1st column: pathology
draw(packLegend(direction = 'vertical', race_lgd, sex_lgd, age_lgd, treatment_lgd), x = unit(85, 'mm'), y = unit(360, 'mm'), just = c('top'))           # 2nd column: demographics
draw(packLegend(direction = 'vertical', blastema_lgd, stroma_lgd, epithelium_lgd, fga_lgd, tmb_lgd), x = unit(135, 'mm'), y = unit(360, 'mm'), just = c('top'))     # 3rd column: genomic events

circos.clear()
graphics.off()

## step 3.8: separate legend-only PDF ####

pdf(file.path(outdir_, 'cohort_summary_circos_legend_only.pdf'), width = 8.5, height = 11)
grid.newpage()
draw(packLegend(direction = 'vertical', gap = unit(2, 'mm'),
                blastema_lgd, stroma_lgd, epithelium_lgd, fga_lgd, tmb_lgd,
                race_lgd, sex_lgd, age_lgd, treatment_lgd,
                survival_lgd, stage_lgd, recurrence_lgd, tumor_lgd, laterality_lgd, size_lgd,
                weight_lgd, histology_lgd, focality_lgd, neoadj_lgd, response_lgd),
     x = unit(0.05, 'npc'), y = unit(0.98, 'npc'), just = c('left', 'top'))
graphics.off()
