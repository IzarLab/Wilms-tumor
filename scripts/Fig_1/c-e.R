# This script plots UMAPs of the whole tumor microenvironment (TME) colored by patient, treatment condition, compartment, subcompartment,
# cell type and cell-cycle phase (one page each)
# Run from the root of the repository: Rscript scripts/Fig_1/c-e.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in TME data ####

load('data/compartments/TME.RData')     # Seurat object tme_

# STEP 2: Plotting TME UMAPs ####

umap_theme = theme(legend.position = 'right',
                   axis.text = element_blank(), axis.ticks = element_blank(),
                   axis.title = element_text(size = 10, face = 'bold'), axis.line = element_line(color = 'black'))
p_ = list()

## step 2.1: patients ####

p_$patient = DimPlot(tme_, group.by = 'orig.ident', raster = F, label = F, shuffle = T, repel = T)+
  umap_theme+
  labs(title = NULL, color = 'Patients', x = 'UMAP', y = 'UMAP')

## step 2.2: treatment condition ####

p_$condition = DimPlot(tme_, group.by = 'condition', raster = F, label = F, repel = T, pt.size = 0.7, shuffle = T, cols = c('red', 'blue'))+
  umap_theme+
  labs(title = NULL, color = 'Treatment', x = 'UMAP', y = 'UMAP')

## step 2.3: compartment ####

p_$compartment = DimPlot(tme_, group.by = 'compartment_name', raster = F, label = F, repel = T, pt.size = 0.7)+
  umap_theme+
  labs(title = NULL, color = 'Compartment', x = 'UMAP', y = 'UMAP')

## step 2.4: subcompartment ####

p_$subcompartment = DimPlot(tme_, group.by = 'subcompartment', raster = F, label = F, repel = F, pt.size = 0.6)+
  umap_theme+
  labs(title = NULL, color = 'Subcompartment', x = 'UMAP', y = 'UMAP')

## step 2.5: cell type (labeled on the plot, no legend) ####

p_$cell_type = DimPlot(tme_, group.by = 'cell_type_2', raster = F, label = T, repel = T, pt.size = 0.7)+
  umap_theme+
  theme(legend.position = 'none', axis.title = element_text(size = 7, face = 'bold'))+
  labs(title = 'Cell type', x = 'UMAP', y = 'UMAP')

## step 2.6: cell-cycle phase ####

p_$phase = DimPlot(tme_, group.by = 'Phase', raster = F, label = F, repel = T)+
  umap_theme+
  labs(title = NULL, color = 'Cell cycle', x = 'UMAP', y = 'UMAP')
rm(tme_); gc()

## step 2.7: saving (one page per UMAP) ####

pdf(file = file.path(outdir_, 'global_umaps.pdf'), width = 8, height = 8)
for(x_ in p_){ plot(x_) }
graphics.off()
rm(p_); gc()
