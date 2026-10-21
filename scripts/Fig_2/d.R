# This script plots the cleaned Jaccard similarity heatmap of NMF factors clustered into programs (metaprograms), annotated by program,
# patient sample and treatment condition
# Input: Jaccard_similarity_matrix_programs.tsv and the manually annotated programs_annotation_processed.tsv (programs_annotation.tsv annotated
# using the program heatmap NMF_programs.pdf and NMF_GSEA_1e-3.tsv)
# NMF terminology: feature = gene; sample = cell; factor = basis component (metagene); rank = max number of factors; V = data matrix (genes x cells);
# W = basis matrix (genes x factors); H = coefficient matrix (factors x cells), each value a metagene expression
# Run from the root of the repository: Rscript scripts/Fig_2/d.R

library(pheatmap)
library(viridis)

set.seed(42)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in data ####

## step 1.1: raw Jaccard similarity matrix ####

mat_ = as.matrix(read.delim(file = 'data/Misc/Jaccard_similarity_matrix_programs.tsv', header = T, sep = '\t', as.is = T, check.names = F))

## step 1.2: manually annotated programs ####

annot_ = read.delim(file = 'data/Misc/programs_annotation_processed.tsv', header = T, sep = '\t', as.is = T, check.names = T)

# STEP 2: Ordering factors by program and condition ####

progs_ = c('blastema 1','S blastema 1','G2M blastema 1',
           'blastema 2',
           'PAX3+ myogenic precursors','striated myocyte',
           'smooth-myocyte precursors','nascent smooth myocyte','smooth myocyte',
           'podocyte',
           'tubular','connecting tubule')     # programs of programs_annotation_processed.tsv, in plotting order
annot_ = do.call(rbind, lapply(progs_, function(p_){ x_ = annot_[annot_$program %in% p_, ]; x_[order(x_$condition), ] }))
rownames(annot_) = annot_$factor_sample
annot_$program = factor(x = annot_$program, levels = progs_)

# STEP 3: Plotting the Jaccard heatmap ####

## step 3.1: annotation colors ####

cols_1 = sample(viridis(n = length(progs_), option = 'turbo'))                    # color for each program
names(cols_1) = levels(annot_$program)
cols_2 = sample(viridis(n = length(unique(annot_$sample)), option = 'turbo'))     # color for each patient sample
names(cols_2) = unique(annot_$sample)
cols_3 = c(untreated = 'blue', neoadjuvant = 'red')                               # color for treatment condition
annot_cols = list(program = cols_1, sample = cols_2, condition = cols_3)

## step 3.2: heatmap (similarities of at most 0.05 set to 0) ####

m_ = mat_[rownames(annot_), rownames(annot_)]
m_[m_ <= 0.05] = 0
rm(mat_)

pdf(file = file.path(outdir_, 'SOM_NMF_programs_cleaned.pdf'), width = 10, height = 10)
pheatmap(mat = m_, annotation_row = annot_[,-1, drop = F], annotation_col = annot_[,-1, drop = F],
         cluster_cols = F, cluster_rows = F, clustering_method = 'ward.D2', treeheight_row = 0, treeheight_col = 0,
         border_color = NA, scale = 'none', color = rev(viridis(n = 100, option = 'A')),
         show_rownames = F, show_colnames = F, fontsize = 2,
         annotation_colors = annot_cols,
         legend_breaks = c(0,.5,1),
         silent = F)
graphics.off()
