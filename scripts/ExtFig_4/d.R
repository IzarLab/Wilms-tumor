# This script plots the plasticity log-ratio, log2((stromal + epithelial + 0.01)/(blastema 1-3 + 0.01)), of malignant cells of treated samples
# with known tumor-burden response by chemotherapy regimen and response (non-responder: unchanged or increased tumor burden; responder:
# decreased)
# Run from the root of the repository: Rscript scripts/ExtFig_4/d.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

set.seed(1)     # reproducible jittered points

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

resp_cols = c('Non-responder' = '#E41A1C', Responder = '#4DAF4A')
resp_labs = c('Non-responder' = 'Non-responder (NR)', Responder = 'Responder (R)')
blastemal_ = c('blastema 1','S blastema 1','G2M blastema 1','blastema 2','G2M blastema 2','blastema 3')
differentiated_ = c('nascent smooth myocyte-like','smooth myocyte-like','S smooth myocyte-like','G2M smooth myocyte-like','striated myocyte-like','myoblast-like',     # stromal
                    'tubules','S tubules','G2M tubules','podocyte-like')                                                                             # epithelial; myogenic precursors are neither

# STEP 1: Malignant cells of treated samples with known tumor-burden response ####

## step 1.1: reading in cancer metadata ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','cell_type')]
rm(s_objs); gc()     # only the metadata is needed

## step 1.2: response group and regimen of each sample (samples of unknown response are dropped) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
cohort_ = cohort_[cohort_$response %in% c('decreased','unchanged','increased'),]
md_ = md_[md_$orig.ident %in% cohort_$patient,]
resp_ = setNames(ifelse(cohort_$response %in% 'decreased', 'Responder', 'Non-responder'), cohort_$patient)
reg_ = setNames(ifelse(cohort_$neoadjuvant %in% 'Vincristine+Doxorubicin', 'V+D', cohort_$neoadjuvant), cohort_$patient)
rm(cohort_)

# STEP 2: Plasticity log-ratio per sample; pseudocount 0.01 avoids inflation when either fraction is near 0 ####

bl_ = tapply(md_$cell_type %in% blastemal_, md_$orig.ident, mean)          # blastemal fraction
dif_ = tapply(md_$cell_type %in% differentiated_, md_$orig.ident, mean)    # differentiated fraction
df_ = data.frame(sample = names(bl_), plasticity = log2((dif_ + 0.01)/(bl_ + 0.01)),
                 response_group = factor(resp_[names(bl_)], levels = names(resp_cols)),
                 regimen = factor(reg_[names(bl_)], levels = c('EE-4A','DD-4A','V+D','I')))
rm(md_, bl_, dif_, resp_, reg_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = regimen, y = plasticity, fill = response_group))+
      theme_classic()+
      theme(axis.text = element_text(color = 'black'), legend.position = 'top')+
      geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8))+
      geom_point(aes(color = response_group), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 2)+
      scale_fill_manual(values = resp_cols, labels = resp_labs, name = 'Response')+
      scale_color_manual(values = resp_cols, labels = resp_labs, name = 'Response')+
      labs(title = 'Plasticity by regimen and response', x = 'Treatment regimen', y = 'Plasticity log-ratio')

ggsave(filename = file.path(outdir_, 'plasticity_log_ratio_by_regimen_and_response.pdf'), plot = p_, width = 7, height = 5)
rm(df_, p_); gc()
