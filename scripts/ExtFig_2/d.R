# This script plots, for each sample, the proportion of each vasculature cell type within the vasculature compartment as a stacked bar plot, with
# cell type colors as in Fig_5/b_d_f.R; chemotherapy-naive samples come first, then chemotherapy-exposed samples
# Run from the root of the repository: Rscript scripts/ExtFig_2/d.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

comp_ = 'vasculature'
comp_cols = c(lymphoid = 'purple3', myeloid = 'orange3', vasculature = 'yellow3')     # compartment colors, as in Fig_5/a.R

# STEP 1: Reading in compartment integrated data ####

load(file.path('data/compartments', paste0('integrated_compartment_', comp_, '_processed.RData')))

# STEP 2: Cell type colors as in Fig_5/b_d_f.R: shades of the compartment color, neighboring cell types (on the UMAP) far apart in shade ####

cen_ = aggregate(Embeddings(s_objs, 'umap')[, 1:2], list(cell_type = as.character(s_objs$cell_type)), median)     # cell type centroids
colnames(cen_)[2:3] = c('x', 'y')
md_ = s_objs@meta.data[, c('orig.ident','condition','cell_type')]
rm(s_objs); gc()     # only the metadata and the centroids are needed

n_ = nrow(cen_)
ord_ = which.min(cen_$x)     # walk from the leftmost cell type to its nearest unvisited neighbor
d_ = as.matrix(dist(cen_[, c('x','y')]))
while(length(ord_) < n_){ r_ = d_[tail(ord_, 1), ]; r_[ord_] = Inf; ord_ = c(ord_, which.min(r_)) }
shade_ = integer(n_)     # along the walk shades alternate between light and dark: 1, n, 2, n-1, ...
shade_[seq(1, n_, 2)] = 1:ceiling(n_/2)
shade_[seq(2, n_, 2)] = n_:(ceiling(n_/2) + 1)
pal_ = colorRampPalette(c('white', comp_cols[[comp_]], 'black'))(n_ + 4)[3:(n_ + 2)]     # from light to dark, without white and black
ct_cols = setNames(pal_[shade_], cen_$cell_type[ord_])
ct_cols = ct_cols[sort(names(ct_cols))]     # alphabetical: first cell type at the bottom of the bars
rm(cen_, n_, ord_, d_, shade_, pal_, r_)

# STEP 3: Proportion of each cell type per sample ####

## step 3.1: percent of cells per sample ####

df_ = as.data.frame(prop.table(table(sample = md_$orig.ident, cell_type = as.character(md_$cell_type)), 1)*100)

## step 3.2: sample order: chemotherapy-naive first, then chemotherapy-exposed ####

cond_ = unique(md_[, c('orig.ident','condition')])
cond_ = cond_[order(factor(cond_$condition, levels = c('untreated','neoadjuvant')), cond_$orig.ident),]
df_$sample = factor(df_$sample, levels = cond_$orig.ident)
df_$cell_type = factor(df_$cell_type, levels = rev(names(ct_cols)))
df_$label = ifelse(df_$Freq >= 5, round(df_$Freq, 1), NA)     # portions below 5% are too thin to be labeled
df_$text_col = ifelse(colSums(col2rgb(ct_cols[as.character(df_$cell_type)])*c(0.299, 0.587, 0.114)) < 128, 'white', 'black')     # white text on dark shades
rm(md_, cond_)

# STEP 4: Plotting ####

p_ =  ggplot(df_, aes(x = sample, y = Freq, fill = cell_type))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'))+
      geom_col(width = 0.8, color = 'black', linewidth = 0.1)+     # thin outline separates neighboring shades
      geom_text(aes(label = label, color = text_col, group = cell_type), position = position_stack(vjust = 0.5), size = 2.3, na.rm = T, show.legend = F)+     # percent at the middle of each portion
      scale_fill_manual(values = ct_cols, breaks = rev(names(ct_cols)), name = NULL)+
      scale_color_identity()+
      scale_y_continuous(expand = c(0,0))+
      labs(title = paste0(tools::toTitleCase(comp_), ' cell types per sample'), x = NULL, y = paste0('% of ', comp_, ' cells'))

ggsave(filename = file.path(outdir_, paste0(comp_, '_cell_type_proportion_per_sample.pdf')), plot = p_, width = 8, height = 4.5)
rm(df_, p_, ct_cols); gc()
