# This script maps malignant cell types to fetal kidney cell types by the overlap between the DGEA markers of each malignant cell type and
# the DGEA markers of each fetal kidney cell type (overlap coefficient and hypergeometric test), as in Fig. 2f
# Run from the root of the repository: Rscript scripts/ExtFig_1/b.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: DGEA between malignant cell types ####

## step 1.1: reading in cancer integrated data; mitotic (S and G2M) cell types are removed ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs[['integrated']] = NULL
s_objs = s_objs[, !grepl('^(G2M|S) ', s_objs$cell_type)]
gc()
genes_cancer = rownames(s_objs[['RNA']])

## step 1.2: markers of each malignant cell type against the others (as in Fig. 2b) ####

Idents(s_objs) = droplevels(as.factor(s_objs$cell_type))
markers_cancer = FindAllMarkers(s_objs, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)
rm(s_objs); gc()
markers_cancer = markers_cancer[markers_cancer$p_val_adj < 0.05,]

# STEP 2: DGEA between fetal kidney cell types ####

## step 2.1: reading in fetal kidney data; mitotic (proliferating) and ureteric bud cell types are removed ####

load('data/compartments/Fetal_kideny_processed.RData')
DefaultAssay(fetal_) = 'RNA'
fetal_ = fetal_[, !grepl('(^Proliferating)|(UB)', fetal_$denovo_cell_type)]
genes_fetal = rownames(fetal_[['RNA']])

## step 2.2: markers of each fetal cell type against the others (as in Fig. 2f) ####

Idents(fetal_) = as.factor(fetal_$denovo_cell_type)
markers_fetal = FindAllMarkers(fetal_, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)
rm(fetal_); gc()
markers_fetal = markers_fetal[markers_fetal$p_val_adj < 0.05,]

# STEP 3: Overlap between malignant and fetal cell type markers ####

## step 3.1: gene universe: genes measured in both cancer and fetal data ####

universe_ = intersect(genes_cancer, genes_fetal)
cancer_markers = lapply(split(markers_cancer$gene, as.character(markers_cancer$cluster)), intersect, universe_)
fetal_markers = lapply(split(markers_fetal$gene, as.character(markers_fetal$cluster)), intersect, universe_)
rm(markers_cancer, markers_fetal, genes_cancer, genes_fetal)

## step 3.2: shared genes, overlap coefficient (shared / size of the smaller set) and hypergeometric test (BH-adjusted) ####

cancer_groups = c('blastema 1' = 'Blastema', 'blastema 2' = 'Blastema', 'blastema 3' = 'Blastema',
                  'smooth myocyte-like precursors' = 'Blastema', 'PAX3+ myogenic precursors' = 'Blastema',
                  'tubules' = 'Epithelium', 'podocyte-like' = 'Epithelium',
                  'nascent smooth myocyte-like' = 'Stroma', 'smooth myocyte-like' = 'Stroma', 'myoblast-like' = 'Stroma', 'striated myocyte-like' = 'Stroma')
ct_order = c('Cap mesenchyme',
             'Proximal renal vesicle','Distal renal vesicle','Proximal S shaped body','Medial S shaped body','Distal S shaped body',
             'Podocyte','Proximal tubule','Loop of Henle',
             'Stroma progenitor','Fibroblast','Mesangial','Myofibroblast','Smooth muscle')
ct_groups = setNames(c('Blastema', rep('Early nephron epithelium', 8), rep('Stroma', 5)), ct_order)
overlap_ = expand.grid(cancer = names(cancer_groups), cell_type = ct_order, stringsAsFactors = F)
overlap_$shared = mapply(function(m_, c_){ length(intersect(cancer_markers[[m_]], fetal_markers[[c_]])) }, overlap_$cancer, overlap_$cell_type)
overlap_$oc = mapply(function(m_, c_, s_){ s_/min(length(cancer_markers[[m_]]), length(fetal_markers[[c_]])) }, overlap_$cancer, overlap_$cell_type, overlap_$shared)
overlap_$p_val = mapply(function(m_, c_, s_){ phyper(s_ - 1, length(fetal_markers[[c_]]), length(universe_) - length(fetal_markers[[c_]]), length(cancer_markers[[m_]]), lower.tail = F) },
                        overlap_$cancer, overlap_$cell_type, overlap_$shared)
overlap_$p_val_adj = p.adjust(overlap_$p_val, method = 'BH')
overlap_$p_val_adj[overlap_$p_val_adj == 0] = .Machine$double.xmin     # p values below ~2e-308 are rounded to 0 (tubules vs Distal S shaped body)
rm(cancer_markers, fetal_markers, universe_)

# STEP 4: Plotting ####

overlap_$cancer = factor(overlap_$cancer, levels = rev(names(cancer_groups)))
overlap_$cancer_group = factor(cancer_groups[as.character(overlap_$cancer)], levels = unique(cancer_groups))
overlap_$cell_type = factor(overlap_$cell_type, levels = ct_order)
overlap_$group = factor(ct_groups[overlap_$cell_type], levels = unique(ct_groups))
overlap_$log_p = -log10(overlap_$p_val_adj)

# setting the color breaks
rng_ = range(overlap_$log_p)
breaks_ = round(seq(from = rng_[1], to = rng_[2], length.out = 4),2)
overlap_$log_p[overlap_$log_p %in% rng_[1]] = breaks_[1]
overlap_$log_p[overlap_$log_p %in% rng_[2]] = breaks_[4]

# plotting
p_ =  ggplot(overlap_, aes(x = cell_type, y = cancer))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_text(size = 7),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1, 'mm'))+
      geom_point(aes(size = oc, color = log_p))+
      scale_color_gradientn(colors = c('grey90','blue','red'), name = '-log10\nadjusted p', breaks = breaks_)+
      scale_radius(range = c(1, 6), name = 'Overlap\ncoefficient')+
      scale_x_discrete(expand = expansion(add = 1))+     # wider panel padding so the one-column Blastema label fits its strip
      facet_grid(cancer_group ~ group, scales = 'free', space = 'free')+
      labs(title = 'Malignant cell types mapped to fetal kidney cell types', x = 'Fetal kidney cell types', y = 'Malignant cell types')

ggsave(filename = file.path(outdir_, 'malignant_cell_types_fetal_kidney_overlap_dotplot.pdf'), plot = p_, width = 8, height = 5.5)
rm(overlap_, p_); gc()
