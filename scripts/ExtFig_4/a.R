# This script plots, for each treated sample with known tumor-burden response, the proportion of malignant cells in each cell group
# (blastema 1-3, myogenic precursors, stromal, epithelial) as a stacked bar plot, split by response (non-responder: unchanged or increased tumor
# burden; responder: decreased tumor burden)
# Run from the root of the repository: Rscript scripts/ExtFig_4/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

grp_cols = c('blastema 1-3' = 'red', 'myogenic precursors' = 'lightgreen', stromal = 'green3', epithelial = 'blue')     # colors of Fig. 2a; precursors color of Fig. 4b
grp_ = c('blastema 1' = 'blastema 1-3', 'S blastema 1' = 'blastema 1-3', 'G2M blastema 1' = 'blastema 1-3', 'blastema 2' = 'blastema 1-3', 'G2M blastema 2' = 'blastema 1-3', 'blastema 3' = 'blastema 1-3',
         'smooth myocyte-like precursors' = 'myogenic precursors', 'G2M smooth myocyte-like precursors' = 'myogenic precursors', 'PAX3+ myogenic precursors' = 'myogenic precursors',
         'nascent smooth myocyte-like' = 'stromal', 'smooth myocyte-like' = 'stromal', 'S smooth myocyte-like' = 'stromal', 'G2M smooth myocyte-like' = 'stromal', 'striated myocyte-like' = 'stromal', 'myoblast-like' = 'stromal',
         'tubules' = 'epithelial', 'S tubules' = 'epithelial', 'G2M tubules' = 'epithelial', 'podocyte-like' = 'epithelial')     # cell group of each cell type

# STEP 1: Malignant cells of treated samples with known tumor-burden response ####

## step 1.1: reading in cancer metadata; cell group of each cell ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','cell_type')]
rm(s_objs); gc()     # only the metadata is needed
md_$group = unname(grp_[as.character(md_$cell_type)])
stopifnot(!anyNA(md_$group))     # every cell type has a group

## step 1.2: response group of each sample (samples of unknown response are dropped) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
cohort_ = cohort_[cohort_$response %in% c('decreased','unchanged','increased'),]
md_ = md_[md_$orig.ident %in% cohort_$patient,]
resp_ = setNames(ifelse(cohort_$response %in% 'decreased', 'Responder', 'Non-responder'), cohort_$patient)
rm(cohort_)

# STEP 2: Proportion of each cell group per sample ####

df_ = as.data.frame(prop.table(table(sample = md_$orig.ident, group = factor(md_$group, levels = names(grp_cols))), 1)*100)
df_$response_group = factor(resp_[as.character(df_$sample)], levels = c('Non-responder','Responder'))

## step 2.1: sample order: by response group, then by decreasing blastemal fraction ####

bl_ = df_[df_$group %in% 'blastema 1-3',]
df_$sample = factor(df_$sample, levels = as.character(bl_$sample[order(bl_$response_group, -bl_$Freq)]))
df_$group = factor(df_$group, levels = rev(names(grp_cols)))     # blastemal at the bottom
df_$label = ifelse(round(df_$Freq, 1) >= 5, round(df_$Freq, 1), NA)     # portions below 5% are too thin to be labeled
df_$text_col = ifelse(df_$group %in% c('blastema 1-3','epithelial'), 'white', 'black')     # white text on red and blue
rm(md_, resp_, bl_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = sample, y = Freq, fill = group))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'))+
      geom_col(color = 'black', linewidth = 0.2)+
      geom_text(aes(label = label, color = text_col, group = group), position = position_stack(vjust = 0.5), size = 2.3, na.rm = T)+     # percent at the middle of each portion
      facet_grid(. ~ response_group, scales = 'free_x', space = 'free_x')+
      scale_fill_manual(values = grp_cols, breaks = names(grp_cols), name = 'Cell group')+
      scale_color_identity()+
      scale_y_continuous(expand = c(0,0))+
      labs(title = 'Malignant composition by tumor-burden response', x = NULL, y = '% of malignant cells')

ggsave(filename = file.path(outdir_, 'cell_group_composition_per_sample_by_response.pdf'), plot = p_, width = 8, height = 5)
rm(df_, p_); gc()
