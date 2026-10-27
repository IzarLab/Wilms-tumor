# This script plots the UMAPs of the fetal kidney nephron and stroma components (kidney cell atlas)
# Run from the root of the repository: Rscript scripts/ExtFig_1/f-g.R

library(Seurat)
options(Seurat.object.assay.version = "v3")
library(ggplot2)
library(ggrepel)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in fetal kidney data ####

## step 1.1: reading in fetal kidney data ####

load('data/compartments/Fetal_kideny_processed.RData')

## step 1.2: modifying cell type annotations; distal UB cell types are removed ####

fetal_ = fetal_[, !fetal_$denovo_cell_type %in% c('CNT/PC - proximal UB','Pelvic epithelium - distal UB')]
fetal_$denovo_cell_type[fetal_$denovo_cell_type %in% 'Proximal UB'] = 'UB'

# STEP 2: Plotting nephron component ####

## step 2.1: data frame of nephron cells and their UMAP coordinates ####

nephron_cells = which(fetal_$component %in% 'nephron')
dt_ = fetal_@meta.data[nephron_cells,]
dt_$x = fetal_[['umap']]@cell.embeddings[nephron_cells,'umap_1']
dt_$y = fetal_[['umap']]@cell.embeddings[nephron_cells,'umap_2']

## step 2.2: cluster centers for labels ####

labels_ = t(sapply(X = unique(dt_$denovo_cell_type), FUN = function(c_){ c(median(dt_$x[dt_$denovo_cell_type %in% c_]), median(dt_$y[dt_$denovo_cell_type %in% c_])) }))
colnames(labels_) = c('x','y')
labels_ = data.frame(labels_, label = rownames(labels_))

## step 2.3: point sizes to make overlapping cell types distinguishable ####

dt_$size = 1
dt_$size[dt_$denovo_cell_type %in% 'Proximal S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Podocyte'] = 4.5

dt_$size[dt_$denovo_cell_type %in% 'Medial S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Proximal tubule'] = 4.5

dt_$size[dt_$denovo_cell_type %in% 'Distal S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Loop of Henle'] = 4.5

## step 2.4: UMAP plot ####

p_ =  ggplot()+
      theme_void()+
      geom_point(data = dt_, aes(x = x, y = y, fill = denovo_cell_type, color = denovo_cell_type), size = dt_$size, show.legend = F, stroke = .2, shape = 21)+
      scale_fill_manual(values = c("Cap mesenchyme" = '#cfd5e8ff',
                                   "Proliferating cap mesenchyme" = '#d3eeeaff',

                                   "Proximal renal vesicle" = '#6b6ca8ff',
                                   "Proximal S shaped body" = '#6b6ca8ff',
                                   "Podocyte" = '#6b6ca8ff',

                                   "Distal renal vesicle" = '#7bb8dcff',

                                   "Medial S shaped body" = '#7bb8dcff',
                                   "Proximal tubule" = '#7bb8dcff',

                                   "Distal S shaped body" = '#96d4d7ff',
                                   "Loop of Henle" = '#96d4d7ff',

                                   "Stroma progenitor" = 'green',
                                   "Proliferating stroma progenitor" = 'purple2',
                                   "Mesangial cells" = 'green4',

                                   "UB" = '#eaa1afff',
                                   "Principal cells" = '#eab1afff',
                                   "Pelvic epithelium" = '#eac1afff'))+
      scale_color_manual(values = c("Cap mesenchyme" = '#cfd5e8ff',
                                    "Proliferating cap mesenchyme" = '#d3eeeaff',

                                    "Proximal renal vesicle" = '#6b6ca8ff',
                                    "Proximal S shaped body" = 'grey20',
                                    "Podocyte" = 'black',

                                    "Distal renal vesicle" = '#7bb8dcff',

                                    "Medial S shaped body" = 'grey40',
                                    "Proximal tubule" = 'black',

                                    "Distal S shaped body" = 'grey40',
                                    "Loop of Henle" = 'black',

                                    "Stroma progenitor" = 'green',
                                    "Proliferating stroma progenitor" = 'purple2',
                                    "Mesangial cells" = 'green4',

                                    "UB" = '#eaa1afff',
                                    "Principal cells" = '#eab1afff',
                                    "Pelvic epithelium" = '#eac1afff'))+
      geom_text_repel(data = labels_, aes(x = x, y = y, label = label), size = 5, color = "black", fontface = 'bold', nudge_y = .6, seed = 3)

ggsave(filename = file.path(outdir_, 'fetal_nephron.pdf'), plot = p_, device = 'pdf', width = 8, height = 8)
rm(dt_, labels_, nephron_cells, p_)

# STEP 3: Plotting stroma component ####

fetal_ = fetal_[, fetal_$component %in% 'stroma']
gc()
p_ = DimPlot(fetal_, group.by = 'denovo_cell_type', pt.size = 1, label = T, label.size = 5, ncol = 1, raster = F, repel = T, shuffle = T) + theme_void() & NoLegend()
ggsave(filename = file.path(outdir_, 'fetal_stroma.pdf'), plot = p_, device = 'pdf', width = 8, height = 8)


