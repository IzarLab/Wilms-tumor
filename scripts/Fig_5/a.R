# This script plots the UMAP of the nonmalignant tumor microenvironment (TME: lymphoid, myeloid and vasculature) colored by subcompartment
# Run from the root of the repository: Rscript scripts/Fig_5/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Subsetting the TME to nonmalignant cells ####

## step 1.1: reading in TME data ####

load('data/compartments/TME.RData')     # Seurat object tme_

## step 1.2: removing malignant cells ####

tme_ = subset(tme_, subset = subcompartment %in% c('lymphoid','myeloid','vasculature'))
tme_$subcompartment = droplevels(factor(tme_$subcompartment))
gc()

# STEP 2: Plotting the nonmalignant TME UMAP ####

## step 2.1: UMAP colored by subcompartment ####

p_ = DimPlot(tme_, group.by = 'subcompartment', raster = F, label = F, repel = T, pt.size = 0.7, cols = c('purple3','orange3','yellow3'))+
  theme(legend.position = 'right',
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))+
  labs(title = NULL, color = 'Compartment', x = 'UMAP', y = 'UMAP')
rm(tme_); gc()

## step 2.2: saving ####

ggsave(filename = file.path(outdir_, 'nonmalignant_TME.pdf'), plot = p_, width = 7, height = 7)
rm(p_); gc()
