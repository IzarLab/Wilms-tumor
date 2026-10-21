# This script finds markers of the malignant subcompartments (blastema, stroma, epithelium) and plots the top 10 markers per subcompartment
# as a dot plot, and the per-sample subcompartment abundance by treatment condition (Wilcoxon test, BH-adjusted)
# Run from the root of the repository: Rscript scripts/Fig_2/b-c.R

library(ggplot2)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(dplyr)
library(tidyr)
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

top10_sub = markers_sub %>%
            filter(!grepl(x = gene, pattern = '(^LINC)|(^A\\w+\\.\\w)')) %>%     # lncRNA exclusion
            filter(pmax(pct.1, pct.2) >= 0.25) %>%
            group_by(cluster) %>%
            arrange(desc(avg_log2FC), .by_group = T) %>%
            slice_head(n = 10) %>%
            ungroup()
write.csv(markers_sub, file.path(outdir_, 'all_markers_subcompartment.csv'), row.names = F)
write.csv(top10_sub, file.path(outdir_, 'top10_markers_per_subcompartment.csv'), row.names = F)
rm(markers_sub)

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

## step 2.4: dot plot (all dots black: dot size alone encodes the percent of cells expressing the gene) ####

p_dot = DotPlot(s_objs, features = gene_order, group.by = 'subcompartment', assay = 'RNA')
p_dot$layers[[1]]$aes_params$colour = 'black'
p_dot = p_dot+
  coord_flip()+
  labs(title = 'Top 10 markers per subcompartment', x = 'Markers', y = 'Subcompartment')+     # genes on the vertical axis after coord_flip()
  guides(size = guide_legend(title = '% of cells\nexpressing'))+
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(color = 'grey90', linewidth = 0.2),
        panel.grid.minor = element_blank())
ggsave(file.path(outdir_, 'dotplot_top10_markers_per_subcompartment_black.pdf'), p_dot, width = 5, height = 8)

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
rm(p_dot, p_box, frac_df, stat_df, ypos_df, top10_sub); gc()
