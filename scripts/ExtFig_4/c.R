# This script plots, for each treated sample with known tumor-burden response, the fraction of malignant cells in each cell type as a
# heatmap, split by response (non-responder: unchanged or increased tumor burden; responder: decreased tumor burden)
# Run from the root of the repository: Rscript scripts/ExtFig_4/c.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

ct_order = c('blastema 1', 'S blastema 1', 'G2M blastema 1',
             'blastema 2', 'G2M blastema 2', 'podocyte-like',                                                     # blastema 2 -> podocyte-like
             'blastema 3', 'smooth myocyte-like precursors', 'G2M smooth myocyte-like precursors', 'nascent smooth myocyte-like',
             'smooth myocyte-like', 'S smooth myocyte-like', 'G2M smooth myocyte-like',                           # blastema 3 -> smooth myocyte-like
             'tubules', 'S tubules', 'G2M tubules',
             'PAX3+ myogenic precursors', 'myoblast-like', 'striated myocyte-like')     # precursors before their progeny (as in ExtFig_3), cycling cells after their cell type

# STEP 1: Malignant cells of treated samples with known tumor-burden response ####

## step 1.1: reading in cancer metadata ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','cell_type')]
rm(s_objs); gc()     # only the metadata is needed

## step 1.2: response group of each sample (samples of unknown response are dropped) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
cohort_ = cohort_[cohort_$response %in% c('decreased','unchanged','increased'),]
md_ = md_[md_$orig.ident %in% cohort_$patient,]
resp_ = setNames(ifelse(cohort_$response %in% 'decreased', 'Responder', 'Non-responder'), cohort_$patient)
rm(cohort_)

# STEP 2: Fraction of each cell type per sample ####

df_ = as.data.frame(prop.table(table(sample = md_$orig.ident, cell_type = factor(md_$cell_type, levels = ct_order)), 1))
df_$response_group = factor(resp_[as.character(df_$sample)], levels = c('Non-responder','Responder'))
df_$cell_type = factor(df_$cell_type, levels = rev(ct_order))     # first cell type at the top

## step 2.1: sample order as in a.R: by response group, then by decreasing blastemal fraction ####

bl_ = tapply(md_$cell_type %in% c('blastema 1','S blastema 1','G2M blastema 1','blastema 2','G2M blastema 2','blastema 3'), md_$orig.ident, mean)
df_$sample = factor(df_$sample, levels = names(bl_)[order(factor(resp_[names(bl_)], levels = c('Non-responder','Responder')), -bl_)])
rm(md_, resp_, bl_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = sample, y = cell_type, fill = Freq))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'))+
      geom_tile(color = 'white', linewidth = 0.2)+
      facet_grid(. ~ response_group, scales = 'free_x', space = 'free_x')+
      scale_fill_gradient(low = 'white', high = 'black', labels = scales::percent_format(), name = 'Fraction')+
      labs(title = 'Malignant cell types by tumor-burden response', x = NULL, y = NULL)

ggsave(filename = file.path(outdir_, 'cell_type_composition_heatmap_by_response.pdf'), plot = p_, width = 8, height = 5.5)
rm(df_, p_); gc()
