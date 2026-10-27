# This scripts plots a heatmap in which programs are matched with annotated cell types using top contributing genes for each program
# NMF-specific terminology:

# Feature: genes
# Sample: cells (in single cell data)
# Factor: basis component (similar to principle component) or a metagene
# Rank: max number of factors (similar to rank in prcomp() for PCA)
# V: data matrix: genes × cells (unlike prcomp that input is cells x genes)
# W: basis matrix (similar to rotation/loadings matrix in PCA): genes × factors
# H: Coefficient Matrix (similar to x matrix in PCA but always non-negative): factors × samples. Unlike PCA, each value is a metagene expression not a coordinate!
# Run from the root of the repository: Rscript scripts/Fig_2/e-f.R

library(dplyr)
library(ggplot2)
library(grid)
library(pheatmap)
library(gridExtra)
library(parallel)
library(doParallel)
library(synchronicity)
library(bigmemory)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(NMF)
library(viridis)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Extracting top contributing genes for each program ####

## step 1.1: reading in table of gene expression programs ####
progs_table = read.delim(file = 'data/Misc/programs_annotation_processed.tsv', header = T, sep = '\t', as.is = T, check.names = T)

## step 1.2: extracting top contributing genes for each program ####
genes_annot = list()
for(m_ in unique(progs_table$program))
{
  if(m_ == 'noise'){ next() }
  message(m_)
  
  factors_ = progs_table$factor_sample[progs_table$program %in% m_]
  genes_ = ws_ = NULL     # genes_ stores top genes and ws_ stores their weights (contributions) for program m
  for(f_ in factors_)
  {
    t_ = strsplit(f_, split = '_|/')[[1]]
    f_num = as.integer(t_[1])     # factor number
    smpl_ = t_[3]                 # sample
    
    # loading NMF results for sample smpl_
    load(paste0('data/Misc/NMF_',smpl_,'.RData'))
    
    # extracting top contributing genes for factor f_num in sample smpl_
    genes_idx = extractFeatures(object = r_estimates$fit[['10']], method = 100)     # 10 is the rank (max number of factors) of the NMF model; method could be a method name or an integer for reporting top genes
    genes_ = c(genes_,featureNames(r_estimates$fit)[genes_idx[[f_num]]])
    
    # extracting weights (contribution) of top genes for factor f_num in sample smpl_
    Ws_ = basis(r_estimates$fit[['10']])                              # basis matrix of sample smpl_
    Ws_ = (Ws_ - min(Ws_))/diff(range(Ws_)) + .Machine$double.eps     # rescaling weights to [0,1] range (NMF does not have a fixed range for weights across different samples)
    ws_fnum = Ws_[ genes_idx[[f_num]], f_num, drop = T]               # weights of all genes for factor f_num in sample smpl_
    ws_ = c(ws_, ws_fnum)
  }
  
  ## step 1.3: filtering top contributing genes for program m_ based on frequency ####
  
  genes_ = sort(table(genes_), decreasing = T)                                                      # frequency of genes across all factors for program m_
  min_freq = if(length(factors_) <= 2){ length(factors_) }else{ ceiling(0.5*length(factors_)) }     # min frequency needed for a top gene to be considered significant
  genes_ = genes_[min_freq <= genes_]                                                               # top genes that are frequent contributors to program m_
  
  ## step 1.4: calculating weights of top genes for program m_ ####
  
  ws_ = sapply(X = names(genes_), FUN = function(g_){ mean(ws_[names(ws_) %in% g_]) })      # weights of top genes for program m_; there could be multiple
                                                                                            # weights for a gene if it is a top contributor to multiple
                                                                                            # factors of program m_
  ## step 1.5: ranking top genes based on frequency and weight for program m_ ####
  
  tmp_ = data.frame(genes = names(genes_), freq = as.numeric(genes_), weights = as.numeric(ws_))
  tmp_ = tmp_[with(tmp_, order(freq, weights, decreasing = T)),]
  tmp_ = tmp_[1:min(5,nrow(tmp_)),]                                                          # top 5 frequent contributing genes to program m_
  
  ## step 1.6: storing top genes and their weights for program m_ ####
  
  genes_annot[[m_]] = data.frame(markers = tmp_$genes, weights = tmp_$weights, frequency = tmp_$freq/length(factors_), program = m_)
}
genes_annot = do.call(genes_annot, what = rbind)

# STEP 2: Computing top gene expressions ####

## step 2.1: reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'

## step 2.2: average expression of each contributing gene in each program across expressing cells of the program's cell type ####

genes_annot$expressions = NA

progs_ = genes_annot$program
progs_[progs_ %in% c("connecting tubule","tubular")] = "tubules"
progs_[progs_ %in% "nascent smooth myocyte"]         = "nascent smooth myocyte-like"
progs_[progs_ %in% "podocyte"]                       = "podocyte-like"
progs_[progs_ %in% "smooth myocyte"]                 = "smooth myocyte-like"
progs_[progs_ %in% "smooth-myocyte precursors"]      = "smooth myocyte-like precursors"
progs_[progs_ %in% "striated myocyte"]               = "striated myocyte-like"

for(m_ in unique(progs_))      # for each program
{
  idx_ = progs_ %in% m_
  genes_ = genes_annot$markers[idx_]
  cells_ = colnames(s_objs)[s_objs$cell_type %in% m_]
  x_ = s_objs[['RNA']]@data[genes_,cells_, drop = F]
  genes_annot$expressions[idx_] = Matrix::rowSums(x_)/pmax(Matrix::rowSums(x_ > 0), 1)     # 0 if no cell expresses the gene
}

# STEP3: Plotting ####

## step 3.1: ordering programs (MP1, MP2, ...) alphabetically and markers by program block ####

prog_levels = sort(unique(genes_annot$program))
prog_levels = prog_levels[order(tolower(prog_levels), method = 'radix')]      # case-insensitive, locale-independent order

plot_df = genes_annot
plot_df$program = factor(plot_df$program, levels = prog_levels)
plot_df = plot_df[order(as.integer(plot_df$program)),]                         # stable: keeps frequency/weight order within each program
plot_df$mp = factor(paste0('MP', as.integer(plot_df$program)), levels = paste0('MP', seq_along(prog_levels)))
plot_df$row_id = factor(paste(plot_df$program, plot_df$markers, sep = '|'),    # a gene can be a marker of several programs
                        levels = rev(unique(paste(plot_df$program, plot_df$markers, sep = '|'))))

## step 3.2: dot plot (gradient = expression, size = NMF weight) ####

# setting the color breaks
rng_ = range(plot_df$expressions)
breaks_ = round(seq(from = rng_[1], to = rng_[2], length.out = 4),2)
plot_df$expressions[plot_df$expressions %in% rng_[1]] = breaks_[1]
plot_df$expressions[plot_df$expressions %in% rng_[2]] = breaks_[4]

# plotting
p_ =  ggplot(plot_df, aes(x = mp, y = row_id, color = expressions, size = weights))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(face = 'italic', size = 8),
            panel.grid.minor = element_blank())+
      geom_point()+
      scale_color_gradientn(colors = c('grey89','blue','red'), name = 'Mean expression in\nexpressing cells', breaks = breaks_)+
      scale_size(range = c(1, 5), name = 'NMF W\nweight')+
      scale_y_discrete(labels = setNames(plot_df$markers, plot_df$row_id))+
      labs(title = 'Top 5 markers per program', x = NULL, y = 'Marker genes ordered by program block')

## step 3.3: program key (MP number -> annotated program) ####

capitalize_ = function(x){ paste0(toupper(substr(x, 1, 1)), substring(x, 2)) }
key_ = textGrob(label = paste0('Annotated programs\n', paste0('MP', seq_along(prog_levels), ': ', capitalize_(prog_levels), collapse = '\n')),
                x = 0, hjust = 0, gp = gpar(fontsize = 10))

## step 3.4: saving ####

ggsave(filename = file.path(outdir_, 'NMF_program_markers_dotplot.pdf'),
       plot = arrangeGrob(p_, key_, ncol = 2, widths = c(3, 1.6)),
       width = 8, height = 10, device = 'pdf')

rm(list = ls()[!ls() %in% c('genes_annot', 'outdir_')]); gc()

# STEP 4: Mapping NMF programs to fetal renal cell types ####

## step 4.1: reading in fetal kidney data ####

load('data/compartments/Fetal_kideny_processed.RData')
fetal_mat = as.matrix(fetal_[['RNA']]@data)
meta_dt = fetal_@meta.data
rm(fetal_); gc()

## step 4.2: calculating average expression and percentage of cells ####

fetal_cell_type = unique(meta_dt$denovo_cell_type)                                          # cell types in fetal kidney data
fetal_cell_type = fetal_cell_type[!grepl(fetal_cell_type, pattern = 'UB')]                  # removing ureteric bud cells from the analysis since they are not part of the nephron lineage.
fetal_cell_type = fetal_cell_type[!grepl(fetal_cell_type, pattern = '^Proliferating')]      # removing mitotic fetal cell types (cell-cycle signal is not lineage-specific)
genes_annot = genes_annot[!grepl(genes_annot$program, pattern = '^(S|G2M) '),]              # removing mitotic NMF programs (S blastema 1, G2M blastema 1)
genes_annot[, paste0(fetal_cell_type,'_expr')] = NA                                         # average expression of each contributing gene in each NMF program across expressing cells of each fetal cell type
genes_annot[, paste0(fetal_cell_type,'_pct')] = NA                                          # percentage of cells expressing each contributing gene in each NMF program across all fetal cell types
for(m_ in unique(genes_annot$program))                                                      # for each NMF program
{
  idx_ = genes_annot$program %in% m_
  genes_ = genes_annot$markers[idx_]
  
  for(ct_ in fetal_cell_type)
  {
    cells_ = rownames(meta_dt)[meta_dt$denovo_cell_type %in% ct_]
    x_ = fetal_mat[genes_, cells_, drop = F]
    genes_annot[idx_, paste0(ct_,'_expr')] = rowSums(x_)/pmax(rowSums(x_ > 0), 1)     # 0 if no cell expresses the gene
    genes_annot[idx_, paste0(ct_,'_pct')] = rowSums(fetal_mat[genes_, cells_, drop = F] > 0)/length(cells_)
  }
}

# STEP 5: Plotting ####

## step 5.1: grouping fetal cell types (columns of the plot) ####

ct_groups = c('Cap mesenchyme'         = 'Blastema',
              'Proximal renal vesicle' = 'Early nephron epithelium',
              'Distal renal vesicle'   = 'Early nephron epithelium',
              'Proximal S shaped body' = 'Early nephron epithelium',
              'Medial S shaped body'   = 'Early nephron epithelium',
              'Distal S shaped body'   = 'Early nephron epithelium',
              'Podocyte'               = 'Early nephron epithelium',
              'Proximal tubule'        = 'Early nephron epithelium',
              'Loop of Henle'          = 'Early nephron epithelium',
              'Stroma progenitor'      = 'Stroma',
              'Fibroblast'             = 'Stroma',
              'Mesangial'              = 'Stroma',
              'Myofibroblast'          = 'Stroma',
              'Smooth muscle'          = 'Stroma')
ct_groups = ct_groups[names(ct_groups) %in% fetal_cell_type]       # only cell types kept in step 4
stopifnot(all(fetal_cell_type %in% names(ct_groups)))               # every fetal cell type must be assigned to a group

## step 5.2: reshaping step 4 results to long format (one row per program, marker and fetal cell type) ####

tmp_ = lapply(names(ct_groups), function(ct_)
                                {
                                  data.frame(program = genes_annot$program,
                                             markers = genes_annot$markers,
                                             cell_type = ct_,
                                             expr = genes_annot[[paste0(ct_,'_expr')]],
                                             pct = genes_annot[[paste0(ct_,'_pct')]])
                                })
genes_dt = do.call(what = rbind, args = tmp_)

## step 5.3: average expression of each program's 5 genes for each fetal cell type and best-matching fetal cell type of each program ####

genes_dt =  genes_dt %>%
            group_by(program, cell_type) %>%
            mutate(average = mean(expr)) %>%          # one value per program and fetal cell type
            group_by(program) %>%
            mutate(best = average == max(average)) %>%
            ungroup() %>%
            as.data.frame()

## step 5.4: ordering rows (programs as in step 3; genes top to bottom followed by the average) and columns ####

progs_ = unique(genes_annot$program)
row_levels = unlist(lapply(progs_, function(m_){ paste(m_, c('Average', rev(genes_annot$markers[genes_annot$program %in% m_])), sep = '|') }))     # bottom to top within each program
genes_dt$row_id = factor(paste(genes_dt$program, genes_dt$markers, sep = '|'), levels = row_levels)
genes_dt$avg_id = factor(paste(genes_dt$program, 'Average', sep = '|'), levels = row_levels)     # row of the program-average square
genes_dt$program = factor(genes_dt$program, levels = progs_)
genes_dt$cell_type = factor(genes_dt$cell_type, levels = names(ct_groups))
genes_dt$group = factor(ct_groups[as.character(genes_dt$cell_type)], levels = unique(ct_groups))

## step 5.5: dot plot (circles: gradient = mean expression in expressing cells, size = % cells; squares: program average, black = best match) ####

genes_dt$pct[genes_dt$pct == 0] = NA     # no circle where no cell expresses the gene
p_ =  ggplot(genes_dt, aes(x = cell_type))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(face = 'italic', size = 8),
            strip.text.x = element_text(size = 7),
            strip.text.y = element_text(angle = 0, hjust = 0),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1, 'mm'))+
      geom_point(aes(y = avg_id), shape = 15, size = 3, color = ifelse(genes_dt$best, 'green2', 'grey70'))+     # program average: same size for all, black = best match; drawn first so the average row stays at the bottom
      geom_point(aes(y = row_id, color = expr, size = 100*pct), na.rm = T)+
      scale_color_gradientn(colors = c('grey89','blue','red'), name = 'Mean expression in\nexpressing cells')+
      scale_radius(range = c(1, 4), limits = c(0, 100), breaks = c(25, 50, 75, 100), name = '% cells')+     # dot radius grows linearly with percent of cells (0% -> 1, 100% -> 4)
      scale_y_discrete(labels = function(x){ sub('^.*\\|', '', x) })+
      facet_grid(program ~ group, scales = 'free', space = 'free', labeller = labeller(group = function(x){ sapply(x, function(v){ paste(strwrap(v, 6), collapse = '\n') }) }))+
      labs(title = 'NMF programs mapped to fetal kidney cell types', x = 'Fetal kidney cell types', y = 'Top 5 NMF program genes')

ggsave(filename = file.path(outdir_, 'NMF_programs_fetal_kidney_dotplot.pdf'), plot = p_, width = 10, height = 13)

