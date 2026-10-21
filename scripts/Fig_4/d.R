# This script builds one Monocle3 trajectory of untreated (chemotherapy-naive) and neoadjuvant-treated (chemotherapy-exposed) malignant cells
# together, rooted at the centroid of the most potent cell type (CytoTRACE 2) as in Fig_3/b-c.R, and shows where along the trajectory
# treated or untreated cells are enriched: the Monocle3 tree is colored by local enrichment and each cell type by its enrichment.
# Plots the tree on the UMAP
# Run from the root of the repository: Rscript scripts/Fig_4/d.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(monocle3)
library(igraph)
library(ggrepel)
# CytoTRACE2 is called as CytoTRACE2::cytotrace2() and not attached: library(CytoTRACE2) also attaches plyr/data.table, which mask dplyr verbs

set.seed(1234)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

enr_cols = c('#2166AC', '#92C5DE', 'grey80', '#F4A582', '#B2182B')      # log2 enrichment of treated over untreated cells: blue = naive, red = exposed
enr_lim = 2                                                             # color scale limits: -enr_lim to enr_lim (4-fold)
pt_cols = c('#2C0E7A', '#5B2A86', '#9A3C8E', '#D85C5C', '#F39B3D', '#F9E721')

# STEP 1: Subsetting the cancer compartment to untreated and treated non-mitotic cells ####

## step 1.1: reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')

## step 1.2: keeping untreated and treated cells and removing mitotic cell types ####

s_objs = s_objs[, s_objs$condition %in% c('untreated','neoadjuvant') & !grepl(s_objs$cell_type, pattern = '^(G2M|S) ')]     # cell-cycle signal is not lineage-specific
gc()
s_objs$cell_type = droplevels(factor(s_objs$cell_type))
print(table(s_objs$cell_type, s_objs$condition))

## step 1.3: sample, condition, cell type and UMAP position of all cells (enrichment uses all cells, not the downsampled ones) ####

all_ = data.frame(sample = s_objs$orig.ident, condition = s_objs$condition, cell_type = as.character(s_objs$cell_type),
                  Embeddings(s_objs, 'umap')[, 1:2])
colnames(all_)[4:5] = c('UMAP_1', 'UMAP_2')

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

## step 2.2: centroid (most central cell) of each cell type and its UMAP position ####

centroid_cells = sapply(core_cells, function(cells_){ cells_[1] })
centroid_xy = Embeddings(s_objs, 'umap')[centroid_cells, 1:2]
rownames(centroid_xy) = names(centroid_cells)

## step 2.3: downsampling (40% of each condition of cell types with at least 2000 cells, all cells otherwise) and dropping all other cells ####

top500_ = sort(unlist(core_cells, use.names = F))
ct_size = table(s_objs$cell_type)     # cell type size, both conditions together: the same fraction is kept from each condition, so treatment ratios do not change
ds_cells = unlist(lapply(split(colnames(s_objs), paste(s_objs$condition, s_objs$cell_type)), function(cells_)
                                                                                             {
                                                                                               if(ct_size[[as.character(s_objs@meta.data[cells_[1], 'cell_type'])]] < 2000){ return(cells_) }     # small cell types: all cells
                                                                                               cells_[sample(length(cells_), round(0.4*length(cells_)))]                         # large cell types: 40%
                                                                                             }), use.names = F)
rm(ct_size)
s_objs = s_objs[, sort(union(top500_, ds_cells))]     # the full object is not needed anymore (frees RAM before CytoTRACE 2)
gc()

# STEP 3: Ordering cell types by developmental potential (CytoTRACE 2) ####

## step 3.1: CytoTRACE 2 on raw counts of the central cells of step 2, in one batch ####

ct_top500 = as.character(s_objs@meta.data[top500_, 'cell_type'])
counts_ = s_objs[['RNA']]@counts[, top500_]
counts_ = counts_[Matrix::rowSums(counts_ > 0) > 0, ]
ct2_ = CytoTRACE2::cytotrace2(as.data.frame(as.matrix(counts_)), species = 'human',
                              batch_size = ncol(counts_), smooth_batch_size = ncol(counts_),     # one batch: the smoothing does not depend on how cells are split
                              parallelize_smoothing = F, ncores = 4)
ct2_ = ct2_[top500_,]
rm(counts_); gc()

## step 3.2: percentage of toti-/pluri-/multipotent cells per cell type ####

potent_ = ct2_$CytoTRACE2_Potency %in% c('Totipotent', 'Pluripotent', 'Multipotent')
pct_potent = sort(tapply(potent_, ct_top500, function(x){ 100*mean(x) }), decreasing = T)
pct_potent = setNames(as.numeric(pct_potent), names(pct_potent))     # named vector, most potent cell type first
print(round(pct_potent, 1))
write.csv(data.frame(cell_type = names(pct_potent), pct_toti_pluri_multipotent = pct_potent),
          file = file.path(outdir_, 'CytoTRACE2_pct_potent_per_celltype_joint.csv'), row.names = F)

# STEP 4: Joint Monocle3 pseudotime of the downsampled cells ####

## step 4.1: keeping the downsampled cells of step 2.3 plus the root cell (centroid of the most potent cell type of step 3) ####

root_ct = names(pct_potent)[1]
root_cell = centroid_cells[[root_ct]]
s_objs = s_objs[, sort(union(ds_cells, root_cell))]
rm(ct2_, top500_, ct_top500, potent_, core_cells, ds_cells); gc()
print(table(s_objs$cell_type, s_objs$condition))

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
cds_ = learn_graph(cds_, use_partition = F, close_loop = F, learn_graph_control = list(minimal_branch_len = 8, ncenter = 800))     # 800 graph nodes so that small cell types get their own branches

## step 4.4: pseudotime rooted at the centroid of the most potent cell type ####

message('Monocle3 root: ', root_cell, ' (', root_ct, ')')
cds_ = order_cells(cds_, root_cells = root_cell)
pt_ = pseudotime(cds_)

df_ = data.frame(cell = colnames(cds_), UMAP_1 = umap_[,1], UMAP_2 = umap_[,2],
                 cell_type = as.character(colData(cds_)$cell_type), pseudotime = as.numeric(pt_[colnames(cds_)]))
rm(umap_, pt_); gc()

# STEP 5: Monocle3 principal graph ####

pg_ = principal_graph(cds_)[['UMAP']]
node_xy = t(cds_@principal_graph_aux[['UMAP']]$dp_mst)[V(pg_)$name, ]     # UMAP position of graph nodes
rm(cds_); gc()

# STEP 6: Enrichment of treated over untreated cells ####

## step 6.1: per cell type: mean share of each cell type across samples of each condition ####

props_ = prop.table(table(all_$sample, all_$cell_type), 1)                 # share of each cell type in each sample
cond_ = tapply(as.character(all_$condition), all_$sample, unique)[rownames(props_)]
ct_enr = data.frame(cell_type = colnames(props_),
                    mean_share_untreated = colMeans(props_[cond_ == 'untreated', , drop = F]),
                    mean_share_treated = colMeans(props_[cond_ == 'neoadjuvant', , drop = F]))
ct_enr$log2_enrichment = log2(ct_enr$mean_share_treated/ct_enr$mean_share_untreated)
rownames(ct_enr) = ct_enr$cell_type
message('samples: ', sum(cond_ == 'untreated'), ' untreated, ', sum(cond_ == 'neoadjuvant'), ' treated')
print(ct_enr)
write.csv(ct_enr, file.path(outdir_, 'joint_trajectory_celltype_enrichment.csv'), row.names = F)
rm(props_, cond_)

## step 6.2: along the tree: share of each condition's cells closest to each graph node, smoothed over graph neighbors ####

all_$node = V(pg_)$name[FNN::get.knnx(node_xy, as.matrix(all_[, c('UMAP_1','UMAP_2')]), k = 1)$nn.index[,1]]     # closest graph node of every cell
all_$w = 1/as.numeric(table(all_$sample)[all_$sample])                                                           # each sample weighs the same
cnt_ = xtabs(w ~ factor(node, levels = V(pg_)$name) + condition, data = all_)
cnt_ = sweep(cnt_, 2, colSums(cnt_), '/')                                                                         # share of each condition
sm_ = t(sapply(ego(pg_, order = 3), function(v_){ colSums(cnt_[as_ids(v_), , drop = F]) }))                       # graph nodes up to 3 edges away
node_enr = setNames(log2((sm_[, 'neoadjuvant'] + 1e-4)/(sm_[, 'untreated'] + 1e-4)), V(pg_)$name)
rm(cnt_, sm_); gc()

# STEP 7: Plotting ####

## step 7.1: tree edges colored by enrichment and cell type centroids ####

el_ = as_edgelist(pg_)
pal_ = colorRampPalette(enr_cols)(101)
tree_ = data.frame(x = node_xy[el_[,1], 1], y = node_xy[el_[,1], 2], xend = node_xy[el_[,2], 1], yend = node_xy[el_[,2], 2],
                   enr = (node_enr[el_[,1]] + node_enr[el_[,2]])/2)
tree_$col = pal_[round((pmin(pmax(tree_$enr, -enr_lim), enr_lim) + enr_lim)/(2*enr_lim)*100) + 1]     # same colors as the fill scale of the centroids

cent_ = data.frame(cell_type = rownames(centroid_xy), x = centroid_xy[,1], y = centroid_xy[,2],
                   n_cells = as.numeric(table(all_$cell_type)[rownames(centroid_xy)]),
                   enr = ct_enr[rownames(centroid_xy), 'log2_enrichment'],
                   label = rownames(centroid_xy), row.names = rownames(centroid_xy))
cent_[root_ct, 'label'] = paste(root_ct, '(root)')

## step 7.2: UMAP with the Monocle3 tree ####

p_umap =  ggplot()+
          geom_point(data = df_[sample(nrow(df_)), ], aes(UMAP_1, UMAP_2, color = pseudotime), size = 0.4, alpha = 0.35, stroke = 0)+
          scale_color_gradientn(colours = pt_cols, name = 'Pseudotime')+
          geom_segment(data = tree_, aes(x = x, y = y, xend = xend, yend = yend), color = tree_$col, linewidth = 1.6, lineend = 'round')+
          geom_point(data = cent_, aes(x, y, size = n_cells, fill = enr), shape = 21, color = 'black', stroke = 1.1)+
          geom_text_repel(data = cent_, aes(x, y, label = label), size = 4.2, point.size = 6, box.padding = 0.5, force = 3, min.segment.length = Inf, seed = 1)+     # labels next to, not on, the circles
          scale_fill_gradientn(colours = enr_cols, limits = c(-enr_lim, enr_lim), oob = scales::squish, breaks = c(-enr_lim, enr_lim), labels = c('naive', 'exposed'), name = 'log2 enrichment')+
          scale_size(range = c(3, 11), name = 'Cells')+
          guides(fill = guide_colorbar(order = 1, title.position = 'top'), size = guide_legend(order = 2, title.position = 'top'),
                 color = guide_colorbar(order = 3, title.position = 'top'))+
          coord_equal()+
          theme_void()+
          theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 18), plot.subtitle = element_text(hjust = 0.5))+
          labs(title = 'Joint trajectory for Chemotherapy-exposed and Chemotherapy-naive',
               subtitle = 'Tree and circles: red = enriched in chemotherapy-exposed, blue = enriched in chemotherapy-naive tumors')

## step 7.3: saving ####

ggsave(filename = file.path(outdir_, 'joint_trajectory_condition_enrichment.pdf'), plot = p_umap, width = 11, height = 11, device = cairo_pdf)
rm(all_, df_, cent_, tree_, p_umap); gc()
