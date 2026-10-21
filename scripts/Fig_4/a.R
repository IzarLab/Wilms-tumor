# This script plots the UMAP of the malignant (cancer) compartment colored by cell type, including mitotic (S/G2M) cell types
# Run from the root of the repository: Rscript scripts/Fig_4/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')

# STEP 2: Plotting the cancer compartment UMAP ####

## step 2.1: cell type colors (mitotic states next to their parent cell type) ####

cell_type_cols = c("blastema 1"                         = 'red4',
                   "S blastema 1"                       = '#460080',
                   "G2M blastema 1"                     = '#C55300',

                   "blastema 2"                         = 'red3',
                   "G2M blastema 2"                     = '#E65300',

                   "blastema 3"                         = 'red',

                   "smooth myocyte-like precursors"     = 'lightgreen',
                   "G2M smooth myocyte-like precursors" = 'green2',

                   "nascent smooth myocyte-like"        = 'green3',

                   "smooth myocyte-like"                = 'green4',
                   "S smooth myocyte-like"              = '#004680',
                   "G2M smooth myocyte-like"            = '#809800',

                   "PAX3+ myogenic precursors"          = '#464600',
                   "myoblast-like"                      = '#E68000',
                   "striated myocyte-like"              = '#DFB25F',

                   "tubules"                            = 'blue',
                   "S tubules"                          = 'lightblue3',
                   "G2M tubules"                        = 'lightblue4',

                   "podocyte-like"                      = 'skyblue')

## step 2.2: UMAP colored and labeled by cell type ####

p_ = DimPlot(s_objs, group.by = 'cell_type', cols = cell_type_cols, pt.size = .4, repel = T, shuffle = T, label = T, raster = F, label.size = 4)+
  theme(legend.position = 'none',
        axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())+
  labs(title = 'Cancer')+
  coord_equal()
q_ = p_; q_$layers[[2]] = NULL     # same UMAP without cell type labels (layer 2 is the label layer)

## step 2.3: saving ####

ggsave(filename = file.path(outdir_, 'malignant_compartment_umap.pdf'), plot = p_, device = 'pdf', width = 7.5, height = 7.5)
ggsave(filename = file.path(outdir_, 'malignant_compartment_umap.png'), plot = q_, device = 'png', width = 7.5, height = 7.5, dpi = 600)
rm(s_objs, p_); gc()
