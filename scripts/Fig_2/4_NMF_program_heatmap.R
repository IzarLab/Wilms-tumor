# This script plots Jaccard similarity matrix of programs (a set of factors defining a cell lineage or state).
# It utilizes GSEA information for a factor to automatically filter out undesired factors. In other words,
# a factor is filtered out, if there is no GSEA associated with it (see NMF_Factors_GSEA.R).

# Input: It needs NMF_metagenes_genes.tsv and NMF_GSEA_1e-3.tsv from NMF_readouts_visualization.R and NMF_Factors_GSEA.R, respectively.
# Output: It outputs a heatmap (NMF_programs.pdf) of clustered factors where each cluster is called a program. Additionally, it outputs
#         programs_annotation.tsv that contains factor-patient pairs in the same order as the clustered rownames (and colnames) in the
#         heatmap. The file programs_annotation.tsv MUST BE ANNOTATED (now named programs_annotation_processed.tsv) USIGING THE OUTPUT
#         HEATMAP (NMF_programs.pdf) AND NMF_GSEA_1e-3.tsv. This way many noisy programs are removed.

# NMF-specific terminology:

# Feature: genes
# Sample: cells (in single cell data)
# Factor: basis component (similar to principle component) or a metagene
# Rank: max number of factors (similar to rank in prcomp() for PCA)
# V: data matrix: genes × cells (unlike prcomp that input is cells x genes)
# W: basis matrix (similar to rotation/loadings matrix in PCA): genes × factors
# H: Coefficient Matrix (similar to x matrix in PCA but always non-negative): factors × samples (cells). Unlike PCA, each value is a metagene expression not a coordinate!

library(pheatmap)
library(gridExtra)
library(viridis)
set.seed(42)

# GSEA of factors per patient ####

ls_ = list()      # for each factor contains top genes
for(f_ in list.files('../out/', pattern = '^W[0-9]+', all.files = T, recursive = T, include.dirs = T))      # for all patient samples
{
  message(f_)
  
  # reading in NMF file
  
  nmf_ = read.delim(file = paste0('../out/',f_,'/NMF/NMF_metagenes_genes.tsv'), header = T, sep = '\t', quote = '', as.is = T, check.names = T)     # top n informative genes defining a factor; padded by a white space because of Excel (from NMF_readouts_visualization.R)
  nmf_gsea = read.delim(file = paste0('../out/',f_,'/NMF/NMF_GSEA_1e-3.tsv'), header = T, sep = '\t', quote = '', as.is = T, check.names = T)       # GSEA on each factor; padded by a white space because of Excel
  
  # GSEA per factor per rank
  # this section excludes factors for which GSEA failed
  
  r_ = 10
  m_ = 100                  # method of reporting top genes
  for(factor_ in 1:r_)      # for each factor (a factor is a metagene)
  {
    cat('\nfactor -> ', factor_,'\n')
    
    gsea_ = nmf_gsea[nmf_gsea$rank %in% r_ & nmf_gsea$factor %in% factor_ & nmf_gsea$nmf_method %in% m_,]     # genes defining this metagene (factor)
    if(nrow(gsea_) == 0){ next() }
    genes_ = trimws(nmf_$genes[nmf_$rank %in% r_ & nmf_$factor %in% factor_ & nmf_$method %in% m_])
    
    ls_[[paste(factor_,f_, sep = '_')]] = genes_
  }
}

# Jaccard similarity matrix ####
# for each pair of factors

mat_ = matrix(data = 0, nrow = length(ls_), ncol = length(ls_), dimnames = list(names(ls_),names(ls_)))
for(r_ in rownames(mat_))
{
  for(c_ in colnames(mat_))
  {
    mat_[r_,c_] = length(intersect(ls_[[r_]], ls_[[c_]]))/length(union(ls_[[r_]], ls_[[c_]]))
  }
}
write.table(x = mat_, file = 'Jaccard_similarity_matrix_programs.tsv', quote = F, sep = '\t', row.names = T, col.names = T)

# Plotting Jaccard heatmap ####

pdf(file = 'NMF_programs.pdf', width = 10, height = 10)
m_ = mat_
m_[m_ <= 0.05] = 0
p_ = pheatmap(mat = m_,
         cluster_cols = T, cluster_rows = T, clustering_method = 'ward.D2', treeheight_row = 0, treeheight_col = 0,
         border_color = NA, scale = 'none', color = rev(viridis(n = 100, option = 'A')),
         show_rownames = T, show_colnames = F, fontsize = 2,
         silent = T)
p_
graphics.off()

# writing annotation file corresponding to clustered heatmap

annot_ = list()
for(prog_ in names(ls_)[p_$tree_row$order])
{
  t_ = strsplit(prog_, split = '_|/', perl = T)[[1]]
  annot_[[prog_]] = data.frame(factor_sample = prog_, program = 'noise', sample = t_[3], condition = t_[2])     # every program is a noise until after manual annotation
}
annot_ = do.call(annot_, what = rbind)
write.table(x = annot_, file = 'programs_annotation.tsv', quote = F, sep = '\t', col.names = T, row.names = F)      # NOTE: ANNOTATE EACH FACTOR USING NMF_programs.pdf AND NMF_GSEA_1e-3.tsv
