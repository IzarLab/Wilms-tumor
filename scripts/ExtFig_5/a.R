# This script plots, for each sample, the proportion of lymphoid, myeloid and vasculature cells within the nonmalignant tumor microenvironment
# (TME) as a stacked bar plot; chemotherapy-naive samples come first, then chemotherapy-exposed samples
# Run from the root of the repository: Rscript scripts/ExtFig_5/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

sub_cols = c(lymphoid = 'purple3', myeloid = 'orange3', vasculature = 'yellow3')     # as in Fig_5/a.R

# STEP 1: Reading in TME metadata ####

load('data/compartments/TME.RData')     # Seurat object tme_
md_ = tme_@meta.data[tme_$subcompartment %in% names(sub_cols), c('orig.ident','condition','subcompartment')]
rm(tme_); gc()     # only the metadata is needed

# STEP 2: Proportion of each subcompartment per sample ####

## step 2.1: percent of cells per sample ####

df_ = as.data.frame(prop.table(table(sample = md_$orig.ident, subcompartment = md_$subcompartment), 1)*100)

## step 2.2: sample order: chemotherapy-naive first, then chemotherapy-exposed ####

cond_ = unique(md_[, c('orig.ident','condition')])
cond_ = cond_[order(factor(cond_$condition, levels = c('untreated','neoadjuvant')), cond_$orig.ident),]
df_$sample = factor(df_$sample, levels = cond_$orig.ident)
df_$subcompartment = factor(df_$subcompartment, levels = rev(names(sub_cols)))     # lymphoid at the bottom
df_$label = ifelse(df_$Freq > 0, round(df_$Freq, 1), NA)     # 0% portions are not labeled
rm(md_, cond_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = sample, y = Freq, fill = subcompartment))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'))+
      geom_col(width = 0.8)+
      geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.3, na.rm = T)+     # percent at the middle of each portion
      scale_fill_manual(values = sub_cols, breaks = names(sub_cols), name = NULL)+
      scale_y_continuous(expand = c(0,0))+
      labs(title = 'Proportion of major TME cell types across all samples', x = NULL, y = 'Proportion within lymphoid + myeloid + vasculature (%)')

ggsave(filename = file.path(outdir_, 'TME_subcompartment_proportion_per_sample.pdf'), plot = p_, width = 8, height = 4.5)
rm(df_, p_); gc()
