# This script plots the UMAPs of the lymphoid, myeloid and vasculature compartments colored and labeled by cell type, with cell types of
# each compartment in shades of the compartment color of a.R
# Run from the root of the repository: Rscript scripts/Fig_5/b_d_f.R

library(ggplot2)
library(Seurat)
library(ggrepel)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

comp_cols = c(lymphoid = 'purple3', myeloid = 'orange3', vasculature = 'yellow3')     # compartment colors, as in a.R

for(comp_ in names(comp_cols))
{
  # STEP 1: Reading in compartment integrated data ####

  load(file.path('data/compartments', paste0('integrated_compartment_', comp_, '_processed.RData')))

  # STEP 2: Plotting the compartment UMAP ####

  ## step 2.1: cell type centroids (median UMAP position) ####

  cen_ = aggregate(Embeddings(s_objs, 'umap')[, 1:2], list(cell_type = as.character(s_objs$cell_type)), median)
  colnames(cen_)[2:3] = c('x', 'y')

  ## step 2.2: cell type colors: shades of the compartment color, neighboring cell types far apart in shade ####

  n_ = nrow(cen_)
  ord_ = which.min(cen_$x)     # walk from the leftmost cell type to its nearest unvisited neighbor
  d_ = as.matrix(dist(cen_[, c('x','y')]))
  while(length(ord_) < n_){ r_ = d_[tail(ord_, 1), ]; r_[ord_] = Inf; ord_ = c(ord_, which.min(r_)) }
  shade_ = integer(n_)     # along the walk shades alternate between light and dark: 1, n, 2, n-1, ...
  shade_[seq(1, n_, 2)] = 1:ceiling(n_/2)
  shade_[seq(2, n_, 2)] = n_:(ceiling(n_/2) + 1)
  pal_ = colorRampPalette(c('white', comp_cols[[comp_]], 'black'))(n_ + 4)[3:(n_ + 2)]     # from light to dark, without white and black
  ct_cols = setNames(pal_[shade_], cen_$cell_type[ord_])
  rm(n_, ord_, d_, shade_, pal_)

  ## step 2.3: UMAP colored by cell type, labels connected to their cell type by a line ####

  p_ = DimPlot(s_objs, group.by = 'cell_type', cols = ct_cols, label = F, pt.size = 1.5)+
    geom_text_repel(data = cen_, aes(x = x, y = y, label = cell_type), inherit.aes = F, size = 5, bg.color = 'white', bg.r = 0.15,
                    box.padding = 1.2, min.segment.length = 0, segment.color = 'black', segment.size = 0.4, max.overlaps = Inf, seed = 1)+
    NoLegend()
  rm(s_objs); gc()

  ## step 2.4: saving ####

  ggsave(filename = file.path(outdir_, paste0(comp_, '_compartment_umap.pdf')), plot = p_, device = 'pdf', width = 7, height = 7)
  rm(p_, cen_, ct_cols); gc()
}
