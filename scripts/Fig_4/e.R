# This script builds one Monocle3 trajectory of neoadjuvant-treated malignant cells of all chemotherapy regimens (EE-4A, DD-4A and others),
# rooted at the centroid of the most potent cell type (CytoTRACE 2) as in Fig_3/b-c.R and d.R, and shows where along the trajectory each
# regimen is enriched: the Monocle3 tree is colored by the dominant regimen and each cell type by the share of each regimen.
# Plots the tree on the UMAP
# Run from the root of the repository: Rscript scripts/Fig_4/e.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(monocle3)
library(igraph)
library(ggrepel)
library(scatterpie)
# CytoTRACE2 is called as CytoTRACE2::cytotrace2() and not attached: library(CytoTRACE2) also attaches plyr/data.table, which mask dplyr verbs

set.seed(1234)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

reg_cols = c(EE4A = '#D73027', DD4A = '#1A9850', Other = '#8C510A')
reg_labs = c(EE4A = 'EE-4A', DD4A = 'DD-4A', Other = 'Other')
pt_cols = c('#2C0E7A', '#5B2A86', '#9A3C8E', '#D85C5C', '#F39B3D', '#F9E721')

# STEP 1: Subsetting the cancer compartment to treated non-mitotic cells ####

## step 1.1: reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')

## step 1.2: keeping treated cells and removing mitotic cell types ####

s_objs = s_objs[, s_objs$condition %in% 'neoadjuvant' & !grepl(s_objs$cell_type, pattern = '^(G2M|S) ')]     # cell-cycle signal is not lineage-specific
gc()
s_objs$cell_type = droplevels(factor(s_objs$cell_type))

## step 1.3: chemotherapy regimen of each sample (EE-4A, DD-4A, and other or unknown regimens) ####

cohort_ = read.delim('data/Misc/cohort.tsv', header = T, sep = '\t', as.is = T)
regimen_ = setNames(ifelse(cohort_$neoadjuvant %in% c('EE-4A','DD-4A'), sub('-', '', cohort_$neoadjuvant), 'Other'), cohort_$patient)
s_objs$regimen = factor(unname(regimen_[as.character(s_objs$orig.ident)]), levels = names(reg_cols))     # unnamed: Seurat matches named vectors to cell names
rm(cohort_, regimen_)
print(table(s_objs$orig.ident, s_objs$regimen))
print(table(s_objs$cell_type, s_objs$regimen))

## step 1.4: sample, regimen, cell type and UMAP position of all cells (used for regimen shares) ####

all_ = data.frame(sample = s_objs$orig.ident, regimen = s_objs$regimen, cell_type = as.character(s_objs$cell_type),
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

## step 2.3: central cells for CytoTRACE 2 (all treated cells are kept for the trajectory: no downsampling) ####

top500_ = sort(unlist(core_cells, use.names = F))

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
          file = file.path(outdir_, 'CytoTRACE2_pct_potent_per_celltype_treated.csv'), row.names = F)

# STEP 4: Monocle3 pseudotime of all treated cells ####

## step 4.1: root cell: centroid of the most potent cell type of step 3 ####

root_ct = names(pct_potent)[1]
root_cell = centroid_cells[[root_ct]]
rm(ct2_, top500_, ct_top500, potent_, core_cells); gc()

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

# STEP 6: Share of each regimen ####

## step 6.1: per cell type: mean share of each cell type across samples of each regimen ####

props_ = prop.table(table(all_$sample, all_$cell_type), 1)                 # share of each cell type in each sample
reg_ = tapply(as.character(all_$regimen), all_$sample, unique)[rownames(props_)]
mean_ = sapply(names(reg_cols), function(r_){ colMeans(props_[reg_ == r_, , drop = F]) })     # cell types x regimens
ct_reg = data.frame(cell_type = rownames(mean_), mean_/rowSums(mean_))                         # regimen shares of each cell type, each regimen weighs the same
message('samples: ', paste(names(table(reg_)), table(reg_), collapse = ', '))
print(ct_reg)
write.csv(data.frame(ct_reg[, 'cell_type', drop = F], setNames(as.data.frame(mean_), paste0('mean_share_', names(reg_cols))), ct_reg[, -1]),
          file.path(outdir_, 'treated_trajectory_celltype_regimen.csv'), row.names = F)
rm(props_, reg_, mean_)

## step 6.2: along the tree: share of each regimen's cells closest to each graph node, smoothed over graph neighbors ####

all_$node = V(pg_)$name[FNN::get.knnx(node_xy, as.matrix(all_[, c('UMAP_1','UMAP_2')]), k = 1)$nn.index[,1]]     # closest graph node of every cell
all_$w = 1/as.numeric(table(all_$sample)[all_$sample])                                                           # each sample weighs the same
cnt_ = xtabs(w ~ factor(node, levels = V(pg_)$name) + regimen, data = all_)
cnt_ = sweep(cnt_, 2, colSums(cnt_), '/')                                                                         # share of each regimen
sm_ = t(sapply(ego(pg_, order = 3), function(v_){ colSums(cnt_[as_ids(v_), , drop = F]) }))                       # graph nodes up to 3 edges away
node_reg = sm_/rowSums(sm_)                                                                                      # regimen shares of each graph node
rownames(node_reg) = V(pg_)$name
rm(cnt_, sm_); gc()

# STEP 7: Plotting ####

## step 7.1: tree edges colored by the dominant regimen and cell type centroids ####

el_ = as_edgelist(pg_)
edge_reg = (node_reg[el_[,1], ] + node_reg[el_[,2], ])/2
tree_ = data.frame(x = node_xy[el_[,1], 1], y = node_xy[el_[,1], 2], xend = node_xy[el_[,2], 1], yend = node_xy[el_[,2], 2],
                   col = reg_cols[max.col(edge_reg)])

r_min = 0.06*diff(range(df_$UMAP_1))*sqrt(min(table(all_$cell_type))/max(table(all_$cell_type)))     # radius of the smallest pie
cent_ = data.frame(cell_type = rownames(centroid_xy), x = centroid_xy[,1], y = centroid_xy[,2],
                   n_cells = as.numeric(table(all_$cell_type)[rownames(centroid_xy)]),
                   ct_reg[rownames(centroid_xy), names(reg_cols)],
                   label = rownames(centroid_xy), row.names = rownames(centroid_xy))
cent_$r = r_min*(cent_$n_cells/min(cent_$n_cells))^0.3     # pie radius grows with the number of treated cells to the power 0.3: large pies cover less of the trajectory
cent_[root_ct, 'label'] = paste(root_ct, '(root)')

leg_ = data.frame(n = c(1000, 5000, 20000))     # pie size legend: pies of 1000, 5000 and 20000 cells side by side at the bottom left
leg_$r = r_min*(leg_$n/min(cent_$n_cells))^0.3
leg_$x = min(df_$UMAP_1) + cumsum(2*leg_$r + 1.2) - leg_$r
leg_$y = min(df_$UMAP_2)

## step 7.2: UMAP with the Monocle3 tree ####

p_umap =  ggplot()+
          geom_point(data = df_[sample(nrow(df_)), ], aes(UMAP_1, UMAP_2, color = pseudotime), size = 0.4, alpha = 0.35, stroke = 0)+
          scale_color_gradientn(colours = pt_cols, name = 'Pseudotime')+
          geom_segment(data = tree_, aes(x = x, y = y, xend = xend, yend = yend), color = tree_$col, linewidth = 1.6, lineend = 'round')+
          geom_scatterpie(data = cent_, aes(x = x, y = y, r = r), cols = names(reg_cols), color = 'black', linewidth = 0.4)+     # pies over the trajectory
          ggforce::geom_circle(data = leg_, aes(x0 = x, y0 = y, r = r), color = 'black', linewidth = 0.4)+     # pie size legend
          geom_text(data = leg_, aes(x = x, y = y - r - 0.3, label = paste(n, 'cells')), size = 3.5)+
          geom_text_repel(data = cent_, aes(x = x, y = y + r, label = label), size = 4.2, nudge_y = 0.25, point.size = NA, min.segment.length = Inf, seed = 1)+
          scale_fill_manual(values = reg_cols, labels = reg_labs, name = 'Regimen')+
          guides(fill = guide_legend(order = 1, title.position = 'top'), color = guide_colorbar(order = 2, title.position = 'top'))+
          coord_equal()+
          theme_void()+
          theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 18), plot.subtitle = element_text(hjust = 0.5))+
          labs(title = 'Joint trajectory for treatment regimens EE-4A vs DD-4A vs Others',
               subtitle = 'Tree: dominant regimen; pies: share of each regimen per cell type, each regimen weighted equally; pie size: number of treated cells')

## step 7.3: saving ####

ggsave(filename = file.path(outdir_, 'treated_trajectory_regimen.pdf'), plot = p_umap, width = 11, height = 11, device = cairo_pdf)
rm(all_, df_, cent_, leg_, tree_, p_umap); gc()
