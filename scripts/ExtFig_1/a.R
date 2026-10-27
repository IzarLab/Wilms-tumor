# This script selects canonical fetal kidney markers from a DGEA between malignant cell types (blastema 1-3 merged, tubules, podocyte-like,
# smooth myocyte-like, striated myocyte-like) and plots their expression in each fetal kidney cell type as a dot plot
# (dot color = mean expression across expressing cells, dot size = percent of cells expressing the gene)
# Run from the root of the repository: Rscript scripts/ExtFig_1/a.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: DGEA between malignant cell types ####

## step 1.1: reading in cancer integrated data and keeping the five (non-mitotic) cell types ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs[['integrated']] = NULL
s_objs = s_objs[, s_objs$cell_type %in% c('blastema 1','blastema 2','blastema 3','tubules','podocyte-like','smooth myocyte-like','striated myocyte-like')]
gc()
s_objs$cell_type = sub('^blastema [1-3]$', 'blastema', as.character(s_objs$cell_type))     # blastema 1-3 merged into blastema

## step 1.2: markers of each cell type against the other four ####

Idents(s_objs) = as.factor(s_objs$cell_type)
markers_ = FindAllMarkers(s_objs, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)
rm(s_objs); gc()

# STEP 2: Canonical fetal markers found by the DGEA ####

## step 2.1: canonical fetal kidney markers (literature) of each group and the malignant cell type each group is tested in ####

canonical_ = list('Blastema-like'                = c('SIX2','CITED1','EYA1','MEOX1','SALL1','PAX2','OSR1','GDNF','CRABP2','UNCX','LEF1'),     # cap mesenchyme / nephron progenitors
                  'Tubular'                      = c('LHX1','PAX8','JAG1','CDH6','HNF1B','EMX2','POU3F3','MECOM','CDH1','KRT8','CUBN','LRP2','SLC34A1','HNF4A','CDH16','SPP1','SLC12A1','UMOD'),     # early nephron and tubules
                  'Podocyte-like'                = c('WT1','TCF21','FOXC2','MAFB','NPHS1','NPHS2','PODXL','PTPRO','CLIC5','PLCE1','SYNPO'),
                  'Smooth myocyte-like'          = c('PDGFRA','PDGFRB','FOXD1','COL1A1','COL3A1','DCN','LUM','POSTN','MYL9','ACTA2','TAGLN','CNN1','MYH11','MYOCD','ACTG2'),
                  'Striated myocyte-like'        = c('MYF5','MYOD1','MYOG','MYMK','DES','CHRNA1','CHRNG','MYH3','MYH8','MYL1','MYL4','ACTC1','ACTA1','TNNT1','TNNT2','TTN','NEB','CKM'))
cell_type_ = c('Blastema-like' = 'blastema','Tubular' = 'tubules', 'Podocyte-like' = 'podocyte-like',
               'Smooth myocyte-like' = 'smooth myocyte-like', 'Striated myocyte-like' = 'striated myocyte-like')

## step 2.2: keeping canonical markers that are significant markers of their malignant cell type, with their highest log2FC in that cell type ####

markers_ = markers_[markers_$p_val_adj < 0.05,]
markers_ = markers_[order(-markers_$avg_log2FC),]
markers_ = markers_[!duplicated(markers_$gene),]     # each gene assigned to the cell type with its highest log2FC
genes_ = lapply(names(canonical_), function(g_)
{
  x_ = markers_[markers_$gene %in% canonical_[[g_]] & markers_$cluster %in% cell_type_[g_],]
  x_$gene[order(match(x_$gene, canonical_[[g_]]))]     # literature order
})
names(genes_) = names(canonical_)
rm(markers_)

# STEP 3: Mean expression and percent of expressing cells per fetal cell type ####

## step 3.1: reading in fetal kidney data; ureteric bud and mitotic (proliferating) cell types are removed ####

load('data/compartments/Fetal_kideny_processed.RData')
genes_ = lapply(genes_, function(g_){ g_[g_ %in% rownames(fetal_[['RNA']]@data)] })     # MYL1, MYMK and MYOD1 are not in the fetal data
mat_ =fetal_[['RNA']]@data[unlist(genes_),]
ct_ = fetal_$denovo_cell_type
rm(fetal_); gc()

blastema_ = c('Cap mesenchyme')
nephron_ = c('Proximal renal vesicle','Distal renal vesicle','Proximal S shaped body','Medial S shaped body',
             'Distal S shaped body','Podocyte','Proximal tubule','Loop of Henle')
stroma_ = c('Stroma progenitor','Fibroblast','Mesangial','Myofibroblast','Smooth muscle')
ct_order = c(blastema_, nephron_, stroma_)     # rows of the plot, top to bottom

## step 3.2: mean of normalized expression (data slot) across expressing cells only, and percent of cells with non-zero expression ####

avg_ = sapply(ct_order, function(c_){ m_ = mat_[, ct_ %in% c_]; Matrix::rowSums(m_)/pmax(Matrix::rowSums(m_ > 0), 1) })     # 0 if no cell expresses the gene
pct_ = sapply(ct_order, function(c_){ Matrix::rowMeans(mat_[, ct_ %in% c_] > 0)*100 })
rm(mat_, ct_)

## step 3.3: long format for plotting ####

df_ = data.frame(gene = rep(rownames(avg_), times = ncol(avg_)),
                 cell_type = rep(colnames(avg_), each = nrow(avg_)),
                 avg_expr = as.vector(avg_), pct = as.vector(pct_))
df_$group = factor(rep(names(genes_), lengths(genes_))[match(df_$gene, unlist(genes_))], levels = names(genes_))
df_$lineage = factor(ifelse(df_$cell_type %in% blastema_, 'Blastema', ifelse(df_$cell_type %in% nephron_, 'Early nephron\nepithelium', 'Stroma')),
                     levels = c('Blastema','Early nephron\nepithelium','Stroma'))
df_$gene = factor(df_$gene, levels = unlist(genes_))
df_$cell_type = factor(df_$cell_type, levels = rev(ct_order))
rm(avg_, pct_)

# STEP 4: Plotting ####

df_$pct[df_$pct == 0] = NA     # no circle where no cell expresses the gene
p_ = ggplot(df_, aes(x = gene, y = cell_type))+
  geom_point(aes(size = pct, color = avg_expr), na.rm = T)+
  facet_grid(lineage ~ group, scales = 'free', space = 'free')+
  scale_color_gradientn(colors = c('grey89','blue','red'), name ='Mean expression in\nexpressing cells')+
  scale_radius(range = c(0, 6), limits = c(0, 100), breaks = c(25, 50, 75), name = '% cells')+     # dot radius proportional to percent of cells, as in Seurat's DotPlot
  labs(title = 'Fetal kidney expression of malignant cell type markers', x = NULL, y = NULL)+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(color = 'black'),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 0.8),
        strip.background = element_rect(fill = 'grey95', color = 'grey60', linewidth = 0.8),
        strip.text = element_text(size = 10),
        strip.text.y = element_text(size = 10, angle = 0),     # horizontal row labels so the one-row blastema panel fits
        legend.position = 'bottom')
ggsave(filename = file.path(outdir_, 'fetal_markers_dotplot.pdf'), plot = p_, width = 18, height = 7)
rm(df_, p_); gc()
