# This script finds markers of the vasculature cell types (DGEA table of all genes) and plots the top 5 protein-coding markers per cell type
# as a dot plot (as in Fig_2/b-c.R)
# Run from the root of the repository: Rscript scripts/ExtFig_2/g.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(dplyr)

if(!requireNamespace('presto', quietly = T)){ stop('presto is required: FindAllMarkers uses its fast Wilcoxon test when installed') }

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

comp_ = 'vasculature'

# STEP 1: Reading in compartment integrated data ####

load(file.path('data/compartments', paste0('integrated_compartment_', comp_, '_processed.RData')))
DefaultAssay(s_objs) = 'RNA'
s_objs[['integrated']] = NULL; gc()     # only the RNA assay is used
s_objs$cell_type = droplevels(factor(s_objs$cell_type))

# STEP 2: Cell type markers ####

## step 2.1: marker discovery on all genes (presto Wilcoxon test) and the full DGEA table ####

Idents(s_objs) = s_objs$cell_type
markers_ = FindAllMarkers(s_objs, only.pos = T, assay = 'RNA', min.pct = 0.09, logfc.threshold = 0.40)
write.csv(markers_, file.path(outdir_, paste0('DGEA_', comp_, '_cell_types.csv')), row.names = F)

## step 2.2: protein-coding genes (for the dot plot only): current and old (alias) HGNC symbols of protein-coding genes in org.Hs.eg.db ####

gt_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGENETYPE)
sym_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
alias_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egALIAS2EG)
pc_ids = gt_$gene_id[gt_$gene_type %in% 'protein-coding']
pc_ = unique(c(sym_$symbol[sym_$gene_id %in% pc_ids], alias_$alias_symbol[alias_$gene_id %in% pc_ids]))
rm(gt_, sym_, alias_, pc_ids)

## step 2.3: top 5 protein-coding markers per cell type, from significant genes detected in at least 25% of cells of the cell type (pct.1) or of all other cells (pct.2) ####

top_ = markers_ %>%
       filter(gene %in% pc_) %>%     # only protein-coding genes (drops lncRNAs and other non-coding genes)
       filter(p_val_adj < 0.05) %>%
       filter(pmax(pct.1, pct.2) >= 0.25) %>%
       group_by(cluster) %>%
       arrange(desc(avg_log2FC), .by_group = T) %>%
       slice_head(n = 5) %>%
       ungroup()
rm(markers_, pc_)

## step 2.4: ordering cell types by hierarchical clustering of their average expression of the top markers ####

avg_mat = AverageExpression(s_objs, assays = 'RNA', group.by = 'cell_type', layer = 'data', return.seurat = F)$RNA
colnames(avg_mat) = levels(s_objs$cell_type)     # AverageExpression replaces some characters of the cell type names
avg_mat = avg_mat[unique(top_$gene), , drop = F]
hc_ = hclust(dist(t(avg_mat)))
group_order = hc_$labels[hc_$order]
rm(avg_mat, hc_)

top_$cluster = factor(top_$cluster, levels = group_order)
gene_order = top_ %>%
             arrange(cluster, desc(avg_log2FC)) %>%
             pull(gene) %>%
             unique() %>%     # a gene that marks two cell types is shown once, with its first cell type
             rev()
rm(top_)

## step 2.5: mean expression in expressing cells and percent of expressing cells per cell type ####

mat_ = s_objs[['RNA']]@data[gene_order,]
ct_ = as.character(s_objs$cell_type)
rm(s_objs); gc()     # the Seurat object is not needed anymore
avg_ = sapply(group_order, function(s_){ x_ = mat_[, ct_ %in% s_]; Matrix::rowSums(x_)/pmax(Matrix::rowSums(x_ > 0), 1) })     # 0 if no cell expresses the gene
pct_ = sapply(group_order, function(s_){ Matrix::rowMeans(mat_[, ct_ %in% s_] > 0)*100 })
rm(mat_, ct_)
dot_df = data.frame(gene = factor(rep(gene_order, times = length(group_order)), levels = gene_order),
                    cell_type = factor(rep(group_order, each = length(gene_order)), levels = group_order),
                    avg_expr = as.vector(avg_), pct = as.vector(pct_))
rm(avg_, pct_)

# STEP 3: Dot plot (dot fill = mean expression in expressing cells, black outline, dot size = percent of cells expressing the gene) ####

dot_df$pct[dot_df$pct == 0] = NA     # no circle where no cell expresses the gene
p_dot = ggplot(dot_df, aes(x = cell_type, y = gene))+
  geom_point(aes(size = pct, fill = avg_expr), shape = 21, color = 'black', na.rm = T)+     # shape 21: filled circle with a separate (black) outline
  scale_fill_gradientn(colors = c('grey89','blue','red'), name = 'Mean expression in\nexpressing cells')+
  scale_radius(range = c(1, 6), limits = c(0, 100), breaks = c(25, 50, 75, 100))+     # dot radius proportional to percent of cells, as in Seurat's DotPlot
  labs(title = paste0('Top 5 markers per ', comp_, ' cell type'), x = 'Cell type', y = 'Markers', size = '% of cells\nexpressing')+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.2),
        panel.grid.minor = element_blank())
ggsave(file.path(outdir_, paste0('dotplot_top5_markers_per_', comp_, '_cell_type.pdf')), p_dot, width = 6, height = 2.5 + 0.13*length(gene_order))     # height grows with the number of genes
rm(dot_df, p_dot); gc()
