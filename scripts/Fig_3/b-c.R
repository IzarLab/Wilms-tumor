# This script orders untreated malignant cell types by developmental potential (CytoTRACE 2), uses the centroid of the most
# potent cell type as root of Monocle3 pseudotime and plots the Monocle3 trajectory tree on the UMAP
# Run from the root of the repository: Rscript scripts/Fig_3/b-c.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(monocle3)
# CytoTRACE2 is called as CytoTRACE2::cytotrace2() and not attached: library(CytoTRACE2) also attaches plyr/data.table, which mask dplyr verbs

set.seed(1234)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Subsetting the cancer compartment to untreated cells ####

## step 1.1: reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')

## step 1.2: keeping untreated cells ####

s_objs = s_objs[, s_objs$condition %in% 'untreated' & !grepl(s_objs$cell_type, pattern = '^(G2M|S) ')]     # removing mitotic cell types (cell-cycle signal is not lineage-specific)
gc()

# STEP 2: Selecting up to 500 most central cells of each cell type ####

## step 2.1: distance of each cell to its cell-type center in PCA space ####

pcs_ = Embeddings(s_objs, 'pca')[, 1:30]     # first 30 PCs
core_cells = lapply(split(colnames(s_objs), as.character(s_objs$cell_type)),  function(cells_)
                                                                              {
                                                                                x_ = pcs_[cells_, , drop = F]
                                                                                d_ = sqrt(rowSums(sweep(x_, 2, apply(x_, 2, median))^2))      # Euclidean distance to the per-PC median of the cell type
                                                                                cells_[order(d_)[1:min(500, length(d_))]]                     # up to 500 closest cells, ordered from the most central
                                                                              })
rm(pcs_); gc()

## step 2.2: centroid (most central cell) of each cell type ####

centroid_cells = sapply(core_cells, function(cells_){ cells_[1] })

## step 2.3: proportional downsampling for step 4 (60% of each cell type with at least 2000 cells, all cells otherwise) and dropping all other cells ####

top500_ = sort(unlist(core_cells, use.names = F))
ds_cells = unlist(lapply(split(colnames(s_objs), as.character(s_objs$cell_type)), function(cells_)
                                                                                  {
                                                                                    if(length(cells_) < 2000){ return(cells_) }                    # small cell types: all cells
                                                                                    cells_[sample(length(cells_), round(0.6*length(cells_)))]      # large cell types: 60% (fills gaps between cell types on the UMAP)
                                                                                  }), use.names = F)
s_objs = s_objs[, sort(union(top500_, ds_cells))]     # the full untreated object is not needed anymore (frees RAM before CytoTRACE 2)
gc()

# STEP 3: Ordering cell types by developmental potential (CytoTRACE 2) ####

## step 3.1: CytoTRACE 2 on raw counts of the central cells of step 2, in one batch ####

ct_top500 = as.character(s_objs@meta.data[top500_, 'cell_type'])
print(table(ct_top500))
counts_ = s_objs[['RNA']]@counts[, top500_]
counts_ = counts_[Matrix::rowSums(counts_ > 0) > 0, ]
ct2_ = CytoTRACE2::cytotrace2(as.data.frame(as.matrix(counts_)), species = 'human',
                              batch_size = ncol(counts_), smooth_batch_size = ncol(counts_),     # one batch: the smoothing does not depend on how cells are split
                              parallelize_smoothing = F, ncores = 4)
ct2_ = ct2_[top500_,]
rm(counts_); gc()

## step 3.2: percentage of toti-/pluri-/multipotent cells per cell type and Kruskal-Wallis test of CytoTRACE 2 scores ####

potent_ = ct2_$CytoTRACE2_Potency %in% c('Totipotent', 'Pluripotent', 'Multipotent')
pct_potent = sort(tapply(potent_, ct_top500, function(x){ 100*mean(x) }), decreasing = T)
pct_potent = setNames(as.numeric(pct_potent), names(pct_potent))     # named vector, most potent cell type first
kw_p = kruskal.test(ct2_$CytoTRACE2_Score, factor(ct_top500))$p.value

print(round(pct_potent, 1))
message('Kruskal-Wallis p-value (CytoTRACE 2 score ~ cell type): ', format.pval(kw_p))
write.csv(data.frame(cell_type = names(pct_potent), pct_toti_pluri_multipotent = pct_potent),
          file = file.path(outdir_, 'CytoTRACE2_pct_potent_per_celltype.csv'), row.names = F)

# STEP 4: Monocle3 pseudotime of the downsampled untreated cells ####

## step 4.1: keeping the downsampled cells of step 2.3 plus the root cell (centroid of the most potent cell type of step 3) ####

root_ct = names(pct_potent)[1]
root_cell = centroid_cells[[root_ct]]
s_objs = s_objs[, sort(union(ds_cells, root_cell))]     # only the root centroid is added: adding other cells changes the small branches of the graph
rm(ct2_, top500_); gc()
print(table(as.character(s_objs$cell_type)))

## step 4.2: cell_data_set from raw counts and the Seurat UMAP ####

counts_ = s_objs[['RNA']]@counts
cds_ = new_cell_data_set(expression_data = counts_, cell_metadata = s_objs@meta.data,
                         gene_metadata = data.frame(gene_short_name = rownames(counts_), row.names = rownames(counts_)))
rm(counts_); gc()
cds_ = cds_[rowSums(counts(cds_) > 0) >= 10, ]     # genes detected in at least 10 cells
umap_ = Embeddings(s_objs, 'umap')[colnames(cds_), 1:2]
colnames(umap_) = c('UMAP_1', 'UMAP_2')
reducedDims(cds_)$UMAP = umap_
rm(s_objs); gc()     # the Seurat object is not needed anymore

## step 4.3: learning the principal graph on the UMAP with cell types as clusters (one partition) ####

cds_@clusters$UMAP$cluster_result$optim_res$membership = as.character(colData(cds_)$cell_type)
cds_@clusters$UMAP$clusters = colData(cds_)$cell_type
cds_@clusters$UMAP$partitions = factor(x = rep(1, ncol(cds_)), levels = '1')
names(cds_@clusters$UMAP$cluster_result$optim_res$membership) = names(cds_@clusters$UMAP$clusters) = names(cds_@clusters$UMAP$partitions) = colnames(cds_)
cds_ = learn_graph(cds_, use_partition = F, close_loop = F, learn_graph_control = list(minimal_branch_len = 10, ncenter = 800))     # 800 graph nodes so that small cell types get their own branches

## step 4.4: pseudotime rooted at the centroid of the most potent cell type of step 3 ####

message('Monocle3 root: ', root_cell, ' (', root_ct, ')')
cds_ = order_cells(cds_, root_cells = root_cell)

# STEP 5: Plotting the Monocle3 trajectory tree on the UMAP ####

root_xy = as.data.frame(reducedDims(cds_)$UMAP[root_cell, , drop = F])     # UMAP position of the root cell
colnames(root_xy) = c('UMAP_1', 'UMAP_2')

## step 5.1: cells colored by pseudotime ####

p_pt = plot_cells(cds_, color_cells_by = 'pseudotime', cell_size = 0.35, trajectory_graph_color = '#00E5FF', trajectory_graph_segment_size = 0.8,     # fluorescent blue: visible on the plasma palette
                  label_cell_groups = F, label_leaves = F, label_branch_points = F, label_roots = F)+
  geom_point(data = root_xy, aes(x = UMAP_1, y = UMAP_2), inherit.aes = F, shape = 21, size = 4, fill = 'red', color = 'black', stroke = 1)+
  annotate('text', x = root_xy$UMAP_1, y = root_xy$UMAP_2, label = 'Root', vjust = -1.2, fontface = 'bold', size = 4)+
  guides(color = guide_colorbar(title = 'Pseudotime'))+
  labs(title = 'Monocle3 trajectory of untreated malignant cells', subtitle = paste0('Root: centroid of ', root_ct, ' (CytoTRACE 2)'))

## step 5.2: cells colored by cell type ####

p_ct = plot_cells(cds_, color_cells_by = 'cell_type', cell_size = 0.35, trajectory_graph_color = 'black', trajectory_graph_segment_size = 0.6,
                  label_cell_groups = F, label_leaves = F, label_branch_points = F, label_roots = F)+
  geom_point(data = root_xy, aes(x = UMAP_1, y = UMAP_2), inherit.aes = F, shape = 21, size = 4, fill = 'red', color = 'black', stroke = 1)+
  annotate('text', x = root_xy$UMAP_1, y = root_xy$UMAP_2, label = 'Root', vjust = -1.2, fontface = 'bold', size = 4)+
  guides(color = guide_legend(title = 'Cell type', override.aes = list(size = 3)))+
  labs(title = 'Monocle3 trajectory of untreated malignant cells', subtitle = paste0('Root: centroid of ', root_ct, ' (CytoTRACE 2)'))

## step 5.3: saving ####

cairo_pdf(filename = file.path(outdir_, 'UMAP_monocle3_trajectory_untreated.pdf'), width = 10, height = 7, onefile = T)
print(p_pt)
print(p_ct)
dev.off()

# STEP 6: Plotting pseudotime by blastema cell type and the other malignant subcompartments ####

## step 6.1: raw Monocle3 pseudotime (same units as the UMAP in step 5.1); blastema split into blastema 1, 2 and 3 ####

pt_ = pseudotime(cds_)
box_dt = data.frame(cell_type = as.character(colData(cds_)$cell_type), subcompartment = as.character(colData(cds_)$subcompartment),
                    pseudotime = pt_)
box_dt$group = ifelse(grepl('^blastema [0-9]$', box_dt$cell_type), box_dt$cell_type, box_dt$subcompartment)
box_dt = box_dt[box_dt$group != 'blastema', ]     # PAX3+ myogenic and smooth myocyte-like precursors (blastema subcompartment, but not blastema 1/2/3) are not shown
box_dt$group = factor(tools::toTitleCase(box_dt$group), levels = c('Blastema 1', 'Blastema 2', 'Blastema 3', 'Epithelium', 'Stroma'))

## step 6.2: boxplot ####

p_box = ggplot(box_dt, aes(x = group, y = pseudotime, fill = group))+
  geom_boxplot(outlier.size = 0.2, width = 0.7, alpha = 0.65, color = 'grey30')+
  scale_fill_manual(values = c('Blastema 1' = '#E64B35', 'Blastema 2' = '#F39B7F', 'Blastema 3' = '#B2331F', Epithelium = '#00A087', Stroma = '#4DBBD5'), guide = 'none')+
  labs(title = 'Pseudotime by malignant subcompartment', x = NULL, y = 'Pseudotime')+
  theme_classic()+
  theme(axis.text = element_text(size = 12, color = 'black'), axis.title.y = element_text(size = 13))

ggsave(filename = file.path(outdir_, 'boxplot_pseudotime_by_subcompartment.pdf'), plot = p_box, width = 7, height = 4.5, device = cairo_pdf)
