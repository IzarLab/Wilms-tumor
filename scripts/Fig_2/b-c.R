# This script finds markers of the malignant subcompartments (blastema, stroma, epithelium) and plots the top 10 markers per subcompartment
# as a dot plot, and the per-sample subcompartment abundance by treatment condition (Wilcoxon test, BH-adjusted)
# Run from the root of the repository: Rscript scripts/Fig_2/b-c.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(dplyr)
library(rstatix)
library(scales)

set.seed(1)     # reproducible jittered points in the box plot

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

condition_cols = c(untreated = '#7F7F7F', neoadjuvant = '#E69F00')

# STEP 1: Reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs$subcompartment = factor(s_objs$subcompartment, levels = c('blastema', 'stroma', 'epithelium'))
s_objs$condition = factor(s_objs$condition, levels = c('untreated', 'neoadjuvant'))

# STEP 2: Subcompartment markers ####

## step 2.1: marker discovery ####

Idents(s_objs) = s_objs$subcompartment
# many canonical blastema genes, including CITED1 and SIX2, are lowly expressed transcription factors that droplet scRNA-seq often fails to
# capture (especially single-nucleus data); they can be biologically real markers and still look like 5-10% of cells in data
markers_sub = FindAllMarkers(s_objs, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)

## step 2.2: top 10 markers per subcompartment, from genes detected in at least 25% of cells of the subcompartment (pct.1) or of all other cells (pct.2) ####

# protein-coding genes: current and old (alias) HGNC symbols of protein-coding genes in org.Hs.eg.db (the data uses some old symbols, e.g. H2AFZ)
gt_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGENETYPE)
sym_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL)
alias_ = AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egALIAS2EG)
pc_ids = gt_$gene_id[gt_$gene_type %in% 'protein-coding']
pc_ = unique(c(sym_$symbol[sym_$gene_id %in% pc_ids], alias_$alias_symbol[alias_$gene_id %in% pc_ids]))
rm(gt_, sym_, alias_, pc_ids)

top10_sub = markers_sub %>%
            filter(gene %in% pc_) %>%     # only protein-coding genes (drops lncRNAs and other non-coding genes)
            filter(pmax(pct.1, pct.2) >= 0.25) %>%
            group_by(cluster) %>%
            arrange(desc(avg_log2FC), .by_group = T) %>%
            slice_head(n = 10) %>%
            ungroup()
write.csv(markers_sub, file.path(outdir_, 'all_markers_subcompartment.csv'), row.names = F)
write.csv(top10_sub, file.path(outdir_, 'top10_markers_per_subcompartment.csv'), row.names = F)
rm(markers_sub, pc_)

## step 2.3: ordering subcompartments by hierarchical clustering of their average expression of the top markers ####

avg_mat = AverageExpression(s_objs, assays = 'RNA', group.by = 'subcompartment', layer = 'data', return.seurat = F)$RNA
avg_mat = avg_mat[intersect(unique(top10_sub$gene), rownames(avg_mat)), , drop = F]
hc_ = hclust(dist(t(avg_mat)))
group_order = hc_$labels[hc_$order]
rm(avg_mat, hc_)

s_objs$subcompartment = factor(s_objs$subcompartment, levels = group_order)
Idents(s_objs) = s_objs$subcompartment
top10_sub$cluster = factor(top10_sub$cluster, levels = group_order)
gene_order = top10_sub %>%
             arrange(cluster, desc(avg_log2FC)) %>%
             distinct(gene, cluster, avg_log2FC) %>%
             pull(gene) %>%
             rev()

## step 2.4: mean expression in expressing cells and percent of expressing cells per subcompartment ####

gene_order = unique(gene_order)
mat_ = s_objs[['RNA']]@data[gene_order,]
avg_ = sapply(group_order, function(s_){ x_ = mat_[, s_objs$subcompartment %in% s_]; Matrix::rowSums(x_)/pmax(Matrix::rowSums(x_ > 0), 1) })     # 0 if no cell expresses the gene
pct_ = sapply(group_order, function(s_){ Matrix::rowMeans(mat_[, s_objs$subcompartment %in% s_] > 0)*100 })
rm(mat_)
dot_df = data.frame(gene = factor(rep(gene_order, times = length(group_order)), levels = gene_order),
                    subcompartment = factor(rep(group_order, each = length(gene_order)), levels = group_order),
                    avg_expr = as.vector(avg_), pct = as.vector(pct_))
rm(avg_, pct_)

## step 2.5: dot plot (dot fill = mean expression in expressing cells, black outline, dot size = percent of cells expressing the gene) ####

dot_df$pct[dot_df$pct == 0] = NA     # no circle where no cell expresses the gene
p_dot = ggplot(dot_df, aes(x = subcompartment, y = gene))+
  geom_point(aes(size = pct, fill = avg_expr), shape = 21, color = 'black', na.rm = T)+     # shape 21: filled circle with a separate (black) outline
  scale_fill_gradientn(colors = c('grey89','blue','red'), name = 'Mean expression in\nexpressing cells')+
  scale_radius(range = c(1, 6), limits = c(0, 100), breaks = c(25, 50, 75, 100))+     # dot radius proportional to percent of cells, as in Seurat's DotPlot
  labs(title = 'Top 10 markers per component', x = 'Component', y = 'Markers', size = '% of cells\nexpressing')+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.2),
        panel.grid.minor = element_blank())
ggsave(file.path(outdir_, 'dotplot_top10_markers_per_subcompartment.pdf'), p_dot, width = 5, height = 8)

# STEP 3: Per-sample subcompartment abundance by condition ####

## step 3.1: fraction of each subcompartment per sample and Wilcoxon test between conditions ####

frac_df = s_objs@meta.data %>%
          select(orig.ident, condition, subcompartment) %>%
          count(orig.ident, condition, subcompartment, name = 'n_cells') %>%
          group_by(orig.ident) %>%
          mutate(total_cells_sample = sum(n_cells), fraction = n_cells/total_cells_sample) %>%
          ungroup()
rm(s_objs); gc()

stat_df = frac_df %>%
          group_by(subcompartment) %>%
          wilcox_test(fraction ~ condition) %>%
          adjust_pvalue(method = 'BH') %>%
          add_significance('p.adj') %>%
          ungroup()
ypos_df = frac_df %>%
          group_by(subcompartment) %>%
          summarise(y.position = max(fraction, na.rm = T)*1.12, .groups = 'drop')
stat_df = stat_df %>%
          left_join(ypos_df, by = 'subcompartment') %>%
          mutate(p_label = paste0(p.adj.signif, ' (BH p = ', signif(p.adj, 3), ')'),
                 xnum = match(subcompartment, levels(frac_df$subcompartment)),
                 xmin = xnum - 0.2, xmax = xnum + 0.2, xmid = xnum)

## step 3.2: box plot with BH-adjusted p-values ####

p_box = ggplot(frac_df, aes(x = subcompartment, y = fraction, fill = condition))+
  geom_boxplot(aes(color = condition), width = 0.7, outlier.shape = NA, alpha = 0.35, position = position_dodge(width = 0.8))+
  geom_point(aes(color = condition), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 2, alpha = 0.85)+
  geom_segment(data = stat_df, aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F, linewidth = 0.4)+
  geom_segment(data = stat_df, aes(x = xmin, xend = xmin, y = y.position*0.985, yend = y.position), inherit.aes = F, linewidth = 0.4)+
  geom_segment(data = stat_df, aes(x = xmax, xend = xmax, y = y.position*0.985, yend = y.position), inherit.aes = F, linewidth = 0.4)+
  geom_text(data = stat_df, aes(x = xmid, y = y.position*1.03, label = p_label), inherit.aes = F, size = 4)+
  scale_fill_manual(values = condition_cols)+
  scale_color_manual(values = condition_cols)+
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.18)))+
  labs(title = 'Subcompartment abundance by condition', x = 'Subcompartment', y = 'Fraction of cells per sample')+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), axis.text.x = element_text(face = 'bold'), legend.title = element_blank())
ggsave(file.path(outdir_, 'boxplot_subcompartment_fraction_by_condition_with_pvalues.pdf'), p_box, width = 8, height = 6)
rm(p_dot, dot_df, p_box, frac_df, stat_df, ypos_df, top10_sub); gc()
