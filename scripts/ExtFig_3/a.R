# This script builds one Monocle3 trajectory of untreated (chemotherapy-naive) and neoadjuvant-treated (chemotherapy-exposed) malignant cells
# together, rooted at the centroid of the most potent cell type (CytoTRACE 2) as in Fig_4/d.R, and plots the scaled pseudotime of each cell
# type by condition as box plots
# Run from the root of the repository: Rscript scripts/ExtFig_3/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(monocle3)
# CytoTRACE2 is called as CytoTRACE2::cytotrace2() and not attached: library(CytoTRACE2) also attaches plyr/data.table, which mask dplyr verbs

set.seed(1234)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

cond_cols = c(untreated = '#2166AC', neoadjuvant = '#B2182B')     # blue = naive, red = exposed (as in Fig_4/d.R)
cond_labs = c(untreated = 'Chemotherapy-naive', neoadjuvant = 'Chemotherapy-exposed')
ct_order = c('blastema 1',
             'blastema 2', 'podocyte-like',                                                                       # blastema 2 -> podocyte-like
             'blastema 3', 'smooth myocyte-like precursors', 'nascent smooth myocyte-like', 'smooth myocyte-like',     # blastema 3 -> smooth myocyte-like
             'tubules',
             'PAX3+ myogenic precursors', 'myoblast-like', 'striated myocyte-like')                               # x-axis: precursors before their progeny

# STEP 1: Subsetting the cancer compartment to untreated and treated non-mitotic cells ####

## step 1.1: reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')

## step 1.2: keeping untreated and treated cells and removing mitotic cell types ####

s_objs = s_objs[, s_objs$condition %in% c('untreated','neoadjuvant') & !grepl(s_objs$cell_type, pattern = '^(G2M|S) ')]     # cell-cycle signal is not lineage-specific
gc()
s_objs$cell_type = droplevels(factor(s_objs$cell_type))
print(table(s_objs$cell_type, s_objs$condition))

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
rm(s_objs, umap_); gc()     # the Seurat object is not needed anymore

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

df_ = data.frame(cell_type = as.character(colData(cds_)$cell_type), condition = as.character(colData(cds_)$condition),
                 pseudotime = as.numeric(pt_[colnames(cds_)]))
rm(cds_, pt_); gc()

# STEP 5: Plotting ####

## step 5.1: pseudotime scaled to [0,1]; cells not reachable from the root (infinite pseudotime) are dropped ####

message('cells with infinite pseudotime: ', sum(!is.finite(df_$pseudotime)))
df_ = df_[is.finite(df_$pseudotime),]
df_$pseudotime = df_$pseudotime/max(df_$pseudotime)
df_$cell_type = factor(df_$cell_type, levels = ct_order)
df_$condition = factor(df_$condition, levels = names(cond_cols))

## step 5.2: box plots of scaled pseudotime by cell type and condition ####

p_ =  ggplot(df_, aes(x = cell_type, y = pseudotime, fill = condition))+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = 'black'),
            axis.text.y = element_text(color = 'black'),
            legend.position = 'right')+
      geom_boxplot(outlier.size = 0.3, linewidth = 0.3)+
      scale_fill_manual(values = cond_cols, labels = cond_labs, name = 'Condition')+
      labs(title = 'Pseudotime by cell type and condition', subtitle = 'Joint trajectory of chemotherapy-naive and chemotherapy-exposed cells', x = NULL, y = 'Scaled pseudotime')

ggsave(filename = file.path(outdir_, 'joint_pseudotime_by_celltype_condition.pdf'), plot = p_, width = 8, height = 5)
rm(df_, p_); gc()
