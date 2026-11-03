# This script plots, for each sample, the proportion of blastema, stroma and epithelium cells within the malignant compartment as a stacked
# bar plot; chemotherapy-naive samples come first, then chemotherapy-exposed samples
# Run from the root of the repository: Rscript scripts/ExtFig_2/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

sub_cols = c(blastema = 'red', stroma = 'green3', epithelium = 'blue')     # as in Fig_2/a.R

# STEP 1: Reading in cancer metadata ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
md_ = s_objs@meta.data[, c('orig.ident','condition','subcompartment')]
rm(s_objs); gc()     # only the metadata is needed

# STEP 2: Proportion of each subcompartment per sample ####

## step 2.1: percent of cells per sample ####

df_ = as.data.frame(prop.table(table(sample = md_$orig.ident, subcompartment = md_$subcompartment), 1)*100)

## step 2.2: sample order: chemotherapy-naive first, then chemotherapy-exposed ####

cond_ = unique(md_[, c('orig.ident','condition')])
cond_ = cond_[order(factor(cond_$condition, levels = c('untreated','neoadjuvant')), cond_$orig.ident),]
df_$sample = factor(df_$sample, levels = cond_$orig.ident)
df_$subcompartment = factor(df_$subcompartment, levels = rev(names(sub_cols)))     # blastema at the bottom
df_$label = ifelse(round(df_$Freq, 1) >= 5, round(df_$Freq, 1), NA)     # portions below 5% are too thin to be labeled
df_$text_col = ifelse(df_$subcompartment %in% 'stroma', 'black', 'white')     # white text on red and blue
rm(md_, cond_)

# STEP 3: Plotting ####

p_ =  ggplot(df_, aes(x = sample, y = Freq, fill = subcompartment))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'),
            legend.position = 'top', legend.justification = 'right')+
      geom_col(width = 0.8)+
      geom_text(aes(label = label, color = text_col, group = subcompartment), position = position_stack(vjust = 0.5), size = 2.3, na.rm = T)+     # percent at the middle of each portion
      scale_fill_manual(values = sub_cols, breaks = names(sub_cols), name = NULL)+
      scale_color_identity()+
      scale_y_continuous(expand = c(0,0))+
      labs(title = 'Malignant components per sample', x = NULL, y = '% of malignant cells')

ggsave(filename = file.path(outdir_, 'malignant_subcompartment_proportion_per_sample.pdf'), plot = p_, width = 8, height = 4.5)
rm(df_, p_); gc()
