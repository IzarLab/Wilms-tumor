# This script computes per-sample plasticity metrics of malignant cells of treated samples with known tumor-burden response, compares
# responders (decreased tumor burden) with non-responders (unchanged or increased tumor burden) by Wilcoxon test and plots the metrics as
# box plots
# Run from the root of the repository: Rscript scripts/ExtFig_4/b.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

set.seed(1)     # reproducible jittered points

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

resp_cols = c('Non-responder' = '#E41A1C', Responder = '#4DAF4A')
grp_ = c('blastema 1' = 'blastemal', 'S blastema 1' = 'blastemal', 'G2M blastema 1' = 'blastemal', 'blastema 2' = 'blastemal', 'G2M blastema 2' = 'blastemal', 'blastema 3' = 'blastemal',
         'smooth myocyte-like precursors' = 'myogenic precursors', 'G2M smooth myocyte-like precursors' = 'myogenic precursors', 'PAX3+ myogenic precursors' = 'myogenic precursors',
         'nascent smooth myocyte-like' = 'stromal', 'smooth myocyte-like' = 'stromal', 'S smooth myocyte-like' = 'stromal', 'G2M smooth myocyte-like' = 'stromal', 'striated myocyte-like' = 'stromal', 'myoblast-like' = 'stromal',
         'tubules' = 'epithelial', 'S tubules' = 'epithelial', 'G2M tubules' = 'epithelial', 'podocyte-like' = 'epithelial')     # cell group of each cell type

# STEP 1: Malignant cells of treated samples with known tumor-burden response ####

## step 1.1: reading in cancer metadata; cell group of each cell ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','cell_type','Phase')]
rm(s_objs); gc()     # only the metadata is needed
md_$group = unname(grp_[as.character(md_$cell_type)])
stopifnot(!anyNA(md_$group))     # every cell type has a group

## step 1.2: response group of each sample (samples of unknown response are dropped) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
cohort_ = cohort_[cohort_$response %in% c('decreased','unchanged','increased'),]
md_ = md_[md_$orig.ident %in% cohort_$patient,]
resp_ = setNames(ifelse(cohort_$response %in% 'decreased', 'Responder', 'Non-responder'), cohort_$patient)
rm(cohort_)

# STEP 2: Plasticity metrics per sample ####

## step 2.1: fraction of each cell group per sample ####

p_ = as.data.frame.matrix(prop.table(table(md_$orig.ident, factor(md_$group, levels = unique(grp_))), 1))

## step 2.2: metrics; differentiated = stromal + epithelial (myogenic precursors count as neither blastemal nor differentiated) ####

met_ = data.frame(sample = rownames(p_), response_group = unname(resp_[rownames(p_)]),
                  blastema = p_$blastemal,
                  differentiated = p_$stromal + p_$epithelial,
                  plasticity = log2((p_$stromal + p_$epithelial + 0.01)/(p_$blastemal + 0.01)),     # pseudocount 0.01 avoids inflation when either fraction is near 0
                  resistance = log2((p_$blastemal + 0.01)/(p_$stromal + p_$epithelial + 0.01)),
                  entropy = -rowSums(p_*log(p_ + 1e-9)))     # Shannon entropy (natural log) of the 4 cell groups
bl_ = md_[md_$group %in% 'blastemal',]
cyc_ = tapply(bl_$Phase %in% c('S','G2M'), bl_$orig.ident, mean)     # fraction of blastemal cells in S or G2M phase
met_$cycling_blastema = ifelse(met_$sample %in% names(cyc_), cyc_[met_$sample], 0)
rm(md_, p_, bl_, cyc_, resp_)

## step 2.3: Wilcoxon test of responders vs non-responders; the unadjusted p-value is used when significant, the BH-adjusted p-value otherwise ####

titles_ = c(blastema = 'Blastema 1-3 persistence', differentiated = 'Total differentiation', plasticity = 'Plasticity log-ratio',
            resistance = 'Resistance log-ratio', entropy = 'State diversity', cycling_blastema = 'Cycling blastema 1-3')
p_val = sapply(names(titles_), function(v_){ wilcox.test(met_[[v_]] ~ met_$response_group, exact = F)$p.value })
p_adj = p.adjust(p_val, method = 'BH')
p_use = ifelse(p_val < 0.05, p_val, p_adj)
sig_ = ifelse(p_use < 0.0001, '****', ifelse(p_use < 0.001, '***', ifelse(p_use < 0.01, '**', ifelse(p_use < 0.05, '*', 'ns'))))
print(data.frame(metric = names(titles_), p_val, p_adj, sig_))

# STEP 3: Plotting ####

## step 3.1: long table of the metrics; panel title = metric (significance) ####

df_ = data.frame(response_group = factor(rep(met_$response_group, length(titles_)), levels = names(resp_cols)),
                 metric = rep(names(titles_), each = nrow(met_)), value = unlist(met_[names(titles_)], use.names = F))
df_$metric = factor(paste0(titles_[df_$metric], ' (', sig_[df_$metric], ')'), levels = paste0(titles_, ' (', sig_, ')'))
rm(met_, titles_, p_val, p_adj, p_use, sig_)

## step 3.2: box plots ####

p_ =  ggplot(df_, aes(x = response_group, y = value, fill = response_group))+
      theme_classic()+
      theme(axis.text = element_text(color = 'black'), legend.position = 'none', strip.background = element_blank())+
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6)+
      geom_jitter(width = 0.12, size = 1.5)+
      facet_wrap(~ metric, nrow = 2, scales = 'free_y')+
      scale_fill_manual(values = resp_cols)+
      scale_x_discrete(labels = c('Non-responder' = 'NR', Responder = 'R'))+
      labs(title = 'Response plasticity metrics', x = NULL, y = NULL)

ggsave(filename = file.path(outdir_, 'plasticity_metrics_by_response.pdf'), plot = p_, width = 8, height = 5)
rm(df_, p_); gc()
