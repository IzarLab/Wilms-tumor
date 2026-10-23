# This script maps malignant subcompartments (blastema, stroma, epithelium) to fetal kidney cell types: for each fetal cell type, its top 6
# DGEA markers (kidney cell atlas, Stewart et al. 2019) are plotted in the malignant subcompartments as a dot plot
# (dot color = mean expression across expressing cells, dot size = percent of cells expressing the gene)
# Run from the root of the repository: Rscript scripts/ExtFig_1/b.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in cancer integrated data and removing mitotic cell types ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs[['integrated']] = NULL
s_objs = s_objs[, !grepl(s_objs$cell_type, pattern = '^(G2M|S) ')]
gc()

# STEP 2: Top 6 markers of each fetal cell type ####

## step 2.1: fetal cell types (rows of the plot, top to bottom); ureteric bud and mitotic (proliferating) cell types are not included ####

nephron_ = c('Cap mesenchyme','Proximal renal vesicle','Distal renal vesicle','Proximal S shaped body','Medial S shaped body',
             'Distal S shaped body','Podocyte','Proximal tubule','Loop of Henle')
stroma_ = c('Stroma progenitor','Fibroblast','Mesangial','Myofibroblast','Smooth muscle')
ct_order = c(nephron_, stroma_)

## step 2.2: top 6 up-regulated markers by log2FC, detected in at least 25% of cells of the fetal cell type and measured in the cancer data ####

fetal_dgea = read.delim(file = 'data/Misc/celltypes_markers_kidneycellatlas.tsv', header = T, sep = '\t', quote = '', as.is = T, check.names = F)
fetal_dgea$gene = trimws(fetal_dgea$gene)
fetal_dgea = fetal_dgea[fetal_dgea$seurat_clusters %in% ct_order & 0 < fetal_dgea$avg_log2FC & fetal_dgea$p_val_adj < 0.05 & 0.25 <= fetal_dgea$pct.1,]
fetal_dgea = fetal_dgea[!grepl(x = fetal_dgea$gene, pattern = '(^LINC)|(^A\\w+\\.\\w)|(-AS\\d)'),]     # lncRNA exclusion
fetal_dgea = fetal_dgea[fetal_dgea$gene %in% rownames(s_objs[['RNA']]@data),]
top6_ = do.call(rbind, lapply(ct_order, function(c_){ x_ = fetal_dgea[fetal_dgea$seurat_clusters %in% c_,]; head(x_[order(-x_$avg_log2FC),], 6) }))
write.csv(top6_, file.path(outdir_, 'fetal_top6_markers_per_cell_type.csv'), row.names = F)
rm(fetal_dgea)

# STEP 3: Mean expression in expressing cells and percent of expressing cells per malignant subcompartment ####

## step 3.1: expression of the markers ####

sub_order = c('blastema','stroma','epithelium')     # columns of the plot
mat_ = s_objs[['RNA']]@data[unique(top6_$gene),]
sub_ = as.character(s_objs$subcompartment)
rm(s_objs); gc()

avg_ = sapply(sub_order, function(s_){ m_ = mat_[, sub_ %in% s_]; Matrix::rowSums(m_)/pmax(Matrix::rowSums(m_ > 0), 1) })     # 0 if no cell expresses the gene
pct_ = sapply(sub_order, function(s_){ Matrix::rowMeans(mat_[, sub_ %in% s_] > 0)*100 })
rm(mat_, sub_)

## step 3.2: long format for plotting (one row per fetal cell type, marker and subcompartment) ####

df_ = do.call(rbind, lapply(sub_order, function(s_)
{
  data.frame(fetal_cell_type = top6_$seurat_clusters, gene = top6_$gene, subcompartment = s_,
             avg_expr = avg_[top6_$gene, s_], pct = pct_[top6_$gene, s_])
}))
df_$fetal_cell_type = factor(df_$fetal_cell_type, levels = ct_order)
df_$subcompartment = factor(df_$subcompartment, levels = sub_order)
df_$gene = factor(df_$gene, levels = rev(unique(top6_$gene)))
write.csv(df_, file.path(outdir_, 'fetal_top6_markers_expression_by_subcompartment.csv'), row.names = F)
rm(avg_, pct_, top6_)

# STEP 4: Plotting ####

p_ = ggplot(df_, aes(x = subcompartment, y = gene))+
  geom_point(aes(size = pct, color = avg_expr))+
  facet_grid(fetal_cell_type ~ ., scales = 'free_y', space = 'free_y')+
  scale_color_gradientn(colors = c('grey89','blue','red'), name = 'Mean expression in\nexpressing cells')+
  scale_size_continuous(range = c(0.5, 5), name = '% cells')+
  labs(title = 'Cancer subsets mapped to fetal kidney', x = NULL, y = NULL)+
  theme_bw(base_size = 11)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black'),
        axis.text.y = element_text(face = 'italic', color = 'black'),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 0.8),
        strip.background = element_rect(fill = 'grey95', color = 'grey60', linewidth = 0.8),
        strip.text.y = element_text(size = 10, angle = 0))
ggsave(filename = file.path(outdir_, 'fetal_markers_by_subcompartment_dotplot.pdf'), plot = p_, width = 6.5, height = 16)
rm(df_, p_); gc()
