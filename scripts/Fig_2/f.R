# This script maps NMF programs to fetal kidney cell types by the overlap between each program's top NMF genes (union of the top 50 genes of
# all its factors) and the DGEA markers of each fetal kidney cell type (overlap coefficient and hypergeometric test)
# Run from the root of the repository: Rscript scripts/Fig_2/f.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(NMF)

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: DGEA between fetal kidney cell types ####

## step 1.1: reading in fetal kidney data; mitotic (proliferating) and ureteric bud cell types are removed ####

load('data/compartments/Fetal_kideny_processed.RData')
DefaultAssay(fetal_) = 'RNA'
fetal_ = fetal_[, !grepl('(^Proliferating)|(UB)', fetal_$denovo_cell_type)]
genes_fetal = rownames(fetal_[['RNA']])

## step 1.2: markers of each fetal cell type against the others (as in Fig. 2b) ####

Idents(fetal_) = as.factor(fetal_$denovo_cell_type)
markers_ = FindAllMarkers(fetal_, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)
rm(fetal_); gc()
markers_ = markers_[markers_$p_val_adj < 0.05,]

# STEP 2: Top NMF genes of each program ####

## step 2.1: factors of non-mitotic programs (noise excluded) ####

progs_table = read.delim(file = 'data/Misc/programs_annotation_processed.tsv', header = T, sep = '\t', as.is = T, check.names = T)
progs_table = progs_table[!progs_table$program %in% 'noise' & !grepl('^(S|G2M) ', progs_table$program),]

## step 2.2: union of the top 50 genes of all factors of each program, and all genes the NMF was run on ####

nmf_genes = list()
nmf_features = NULL
for(f_ in list.files('data/Misc', pattern = '^NMF_.*\\.RData$', full.names = T))     # one NMF result per sample
{
  load(f_)
  smpl_ = sub('^NMF_(.*)\\.RData$', '\\1', basename(f_))
  nmf_features = union(nmf_features, featureNames(r_estimates$fit))
  top_ = extractFeatures(object = r_estimates$fit[['10']], method = 50)     # top 50 genes of each factor (rank 10 model)
  for(i_ in which(sub('^.*/', '', progs_table$factor_sample) %in% smpl_))     # factors of this sample assigned to a program
  {
    f_num = as.integer(strsplit(progs_table$factor_sample[i_], split = '_')[[1]][1])
    m_ = progs_table$program[i_]
    nmf_genes[[m_]] = union(nmf_genes[[m_]], featureNames(r_estimates$fit)[top_[[f_num]]])
  }
  rm(r_estimates)
}

# STEP 3: Overlap between program genes and fetal cell type markers ####

## step 3.1: gene universe: genes the NMF was run on that are also measured in fetal data ####

universe_ = intersect(nmf_features, genes_fetal)
nmf_genes = lapply(nmf_genes, intersect, universe_)
fetal_markers = lapply(split(markers_$gene, as.character(markers_$cluster)), intersect, universe_)
rm(markers_, nmf_features, genes_fetal)

## step 3.2: shared genes, overlap coefficient (shared / size of the smaller set) and hypergeometric test (BH-adjusted) ####

progs_ = unique(progs_table$program)     # program order as in the annotation table
ct_order = c('Cap mesenchyme',
             'Proximal renal vesicle','Distal renal vesicle','Proximal S shaped body','Medial S shaped body','Distal S shaped body',
             'Podocyte','Proximal tubule','Loop of Henle',
             'Stroma progenitor','Fibroblast','Mesangial','Myofibroblast','Smooth muscle')
ct_groups = setNames(c('Blastema', rep('Early nephron epithelium', 8), rep('Stroma', 5)), ct_order)
overlap_ = expand.grid(program = progs_, cell_type = ct_order, stringsAsFactors = F)
overlap_$shared = mapply(function(p_, c_){ length(intersect(nmf_genes[[p_]], fetal_markers[[c_]])) }, overlap_$program, overlap_$cell_type)
overlap_$oc = mapply(function(p_, c_, s_){ s_/min(length(nmf_genes[[p_]]), length(fetal_markers[[c_]])) }, overlap_$program, overlap_$cell_type, overlap_$shared)
overlap_$p_val = mapply(function(p_, c_, s_){ phyper(s_ - 1, length(fetal_markers[[c_]]), length(universe_) - length(fetal_markers[[c_]]), length(nmf_genes[[p_]]), lower.tail = F) },
                        overlap_$program, overlap_$cell_type, overlap_$shared)
overlap_$p_val_adj = p.adjust(overlap_$p_val, method = 'BH')
rm(nmf_genes, fetal_markers, universe_)

# STEP 4: Plotting ####

overlap_$program = factor(overlap_$program, levels = rev(progs_))
overlap_$cell_type = factor(overlap_$cell_type, levels = ct_order)
overlap_$group = factor(ct_groups[overlap_$cell_type], levels = unique(ct_groups))
overlap_$log_p = -log10(overlap_$p_val_adj)

# setting the color breaks
rng_ = range(overlap_$log_p)
breaks_ = round(seq(from = rng_[1], to = rng_[2], length.out = 4),2)
overlap_$log_p[overlap_$log_p %in% rng_[1]] = breaks_[1]
overlap_$log_p[overlap_$log_p %in% rng_[2]] = breaks_[4]

# plotting
p_ =  ggplot(overlap_, aes(x = cell_type, y = program))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.x = element_text(size = 7),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1, 'mm'))+
      geom_point(aes(size = oc, color = log_p))+
      scale_color_gradientn(colors = c('grey90','blue','red'), name = '-log10\nadjusted p', breaks = breaks_)+
      scale_radius(range = c(1, 6), name = 'Overlap\ncoefficient')+
      scale_x_discrete(expand = expansion(add = 1))+     # wider panel padding so the one-column Blastema label fits its strip
      facet_grid(. ~ group, scales = 'free_x', space = 'free_x')+
      labs(title = 'NMF programs mapped to fetal kidney cell types', x = 'Fetal kidney cell types', y = 'NMF programs')

ggsave(filename = file.path(outdir_, 'NMF_programs_fetal_kidney_overlap_dotplot.pdf'), plot = p_, width = 8, height = 5)
