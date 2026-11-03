# This script plots the blastema 1-3 fraction (blastema 1-3 and their mitotic cells) of malignant cells of treated samples with known tumor-burden
# response by response category (decreased, unchanged and increased tumor burden)
# Run from the root of the repository: Rscript scripts/ExtFig_4/f.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

set.seed(1)     # reproducible jittered points

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

resp_labs = c(decreased = 'Decreased', unchanged = 'Unchanged', increased = 'Increased')
blastemal_ = c('blastema 1','S blastema 1','G2M blastema 1','blastema 2','G2M blastema 2','blastema 3')

# STEP 1: Malignant cells of treated samples with known tumor-burden response ####

## step 1.1: reading in cancer metadata ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','cell_type')]
rm(s_objs); gc()     # only the metadata is needed

## step 1.2: response category of each sample (samples of unknown response are dropped) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
cohort_ = cohort_[cohort_$response %in% names(resp_labs),]
md_ = md_[md_$orig.ident %in% cohort_$patient,]
resp_ = setNames(cohort_$response, cohort_$patient)
rm(cohort_)

# STEP 2: Blastemal fraction per sample ####

bl_ = tapply(md_$cell_type %in% blastemal_, md_$orig.ident, mean)
df_ = data.frame(sample = names(bl_), blastema = as.numeric(bl_), response = factor(resp_[names(bl_)], levels = names(resp_labs)))
rm(md_, bl_, resp_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = response, y = blastema, fill = response))+
      theme_classic()+
      theme(axis.text = element_text(color = 'black'), legend.position = 'none')+
      geom_boxplot(outlier.shape = NA, alpha = 0.7)+
      geom_jitter(width = 0.12, size = 2)+
      scale_x_discrete(labels = resp_labs)+
      labs(title = 'Blastema 1-3 persistence by tumor burden', x = 'Tumor burden', y = 'Blastema 1-3 fraction')

ggsave(filename = file.path(outdir_, 'blastema_fraction_by_response_category.pdf'), plot = p_, width = 5.5, height = 4.5)
rm(df_, p_); gc()
