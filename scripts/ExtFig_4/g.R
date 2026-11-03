# This script plots the mean, over treated samples, of the proportion of malignant cells in each cell group (blastema 1-3, myogenic precursors,
# stromal, epithelial) for non-responders (unchanged or increased tumor burden) and responders (decreased tumor burden) as stacked bars
# Run from the root of the repository: Rscript scripts/ExtFig_4/g.R

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

# STEP 2: Mean proportion of each cell group per response group ####

p_ = prop.table(table(md_$orig.ident, factor(md_$group, levels = names(grp_cols))), 1)*100     # per sample, 0 where a sample has no cell of a group
mean_ = apply(p_, 2, function(x){ tapply(x, resp_[rownames(p_)], mean) })     # response group x cell group
df_ = data.frame(response_group = factor(rep(rownames(mean_), ncol(mean_)), levels = c('Non-responder','Responder')),
                 group = factor(rep(colnames(mean_), each = nrow(mean_)), levels = rev(names(grp_cols))),     # blastemal at the bottom
                 mean_prop = as.vector(mean_))
df_$text_col = ifelse(df_$group %in% c('blastema 1-3','epithelial'), 'white', 'black')     # white text on red and blue
print(df_)
rm(md_, resp_, p_, mean_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = response_group, y = mean_prop, fill = group))+
      theme_classic()+
      theme(axis.text = element_text(color = 'black'))+
      geom_col(color = 'black', linewidth = 0.25, width = 0.75)+
      geom_text(aes(label = round(mean_prop, 1), color = text_col, group = group), position = position_stack(vjust = 0.5), size = 3)+     # mean percent at the middle of each portion
      scale_fill_manual(values = grp_cols, breaks = names(grp_cols), name = 'Cell group')+
      scale_color_identity()+
      scale_x_discrete(labels = c('Non-responder' = 'NR', Responder = 'R'))+
      scale_y_continuous(expand = c(0,0))+
      labs(title = 'Mean composition by response', x = NULL, y = 'Mean % of malignant cells')

ggsave(filename = file.path(outdir_, 'mean_cell_group_composition_by_response.pdf'), plot = p_, width = 4.5, height = 5)
rm(df_, p_); gc()
