# This script plots the UMAP of the malignant (cancer) compartment colored by subcompartment (blastema, stroma, epithelium)
# Run from the root of the repository: Rscript scripts/Fig_2/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs$subcompartment = factor(s_objs$subcompartment, levels = c('blastema', 'stroma', 'epithelium'))

# STEP 2: Plotting the subcompartment UMAP ####

## step 2.1: UMAP colored by subcompartment ####

p_ = DimPlot(s_objs, group.by = 'subcompartment', cols = c(blastema = 'red', epithelium = 'blue', stroma = 'green3'), raster = F, label = F, repel = F, pt.size = 0.6)+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Subcompartment', x = 'UMAP', y = 'UMAP')
rm(s_objs); gc()

## step 2.2: saving ####

ggsave(filename = file.path(outdir_, 'WT_component_umap.pdf'), plot = p_, width = 8, height = 8)
rm(p_); gc()
