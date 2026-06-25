# ============================================================
# Figure 2 fetal kidney reference marker dot plots
# ============================================================
#
# GitHub script name : 03_plot_fetal_kidney_reference_markers.R
# Original file      : Fig.2.fetal_nephron.R
# Manuscript context : Wilms tumor developmental plasticity project
#
# Purpose:
#   Plots fetal kidney/nephron and stromal reference marker expression used for developmental-state interpretation.
#
# Usage:
#   Run from the root of the GitHub repository, or edit the input/output
#   paths in the SETTINGS / Input-output section below before running.
#
# Notes for reproducibility:
#   - This script preserves the analytical logic of the original working
#     manuscript script.
#   - Hard-coded local paths from the original analysis may need to be
#     changed to repository-relative paths such as data/processed/ and
#     results/figure_XX/.
#   - Large objects such as Seurat .RData files should usually be stored
#     outside GitHub or tracked with Git LFS / external data release.
#
# Suggested output folder:
#   results/figure_02/
#
# ============================================================

# This scripts plots fetal nephron and stroma from kidney cell atlas

library(Seurat)
#options(Seurat.object.assay.version = "v3")
library(SeuratObject)
#library(SeuratWrappers)
library(ggplot2)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(cowplot)

# Reading compartments data ####

obj<-load(file = '~/Documents/Izar_Group/wilms/Fetal_kideny_processed.RData')
obj
## modifying cell type annotations ####

# same cleanup as before
fetal_ <- subset(
  fetal_,
  cells = colnames(fetal_)[!fetal_$denovo_cell_type %in% c(
    "CNT/PC - proximal UB",
    "Pelvic epithelium - distal UB"
  )]
)
fetal_$denovo_cell_type[fetal_$denovo_cell_type %in% "Proximal UB"] <- "UB"

DefaultAssay(fetal_) <- "RNA"
Idents(fetal_) <- "denovo_cell_type"

# order the x-axis groups explicitly
group_order <- c(
  "Cap mesenchyme",
  "Proliferating cap mesenchyme",
  "Proximal renal vesicle",
  "Distal renal vesicle",
  "Proximal S shaped body",
  "Medial S shaped body",
  "Distal S shaped body",
  "Podocyte",
  "Proximal tubule",
  "Loop of Henle",
  "UB",
  "Stroma progenitor",
  "Proliferating stroma progenitor",
  "Fibroblast",
  "Proliferating fibroblast",
  "Mesangial",
  "Myofibroblast",
  "Smooth muscle"
)

group_order <- group_order[group_order %in% unique(fetal_$denovo_cell_type)]
fetal_$denovo_cell_type <- factor(fetal_$denovo_cell_type, levels = group_order)
Idents(fetal_) <- "denovo_cell_type"

# markers in ordered blocks
marker_blocks <- list(
  "Cap mesenchyme / blastema-like" = c("SIX2", "CITED1", "EYA1", "WT1", "SALL1", "PAX2"),
  "Early nephron epithelium"       = c("PAX8", "EPCAM", "CDH6"),
  "Tubular"                        = c("LRP2", "SLC12A1"),
  "Podocyte-like"                  = c("NPHS1", "PODXL"),
  "Stromal / smooth muscle"        = c("TAGLN", "MYL9", "PDGFRA"),
  "Myogenic"                       = c("MYOD1", "MYOG", "MYH3")
)

# keep only genes present in object
marker_blocks <- lapply(marker_blocks, function(x) x[x %in% rownames(fetal_)])
marker_blocks <- marker_blocks[sapply(marker_blocks, length) > 0]

# single dotplot
p <- DotPlot(
  object = fetal_,
  features = marker_blocks,
  group.by = "denovo_cell_type",
  cols = c("lightgrey", "black")
) +
  RotatedAxis() +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text.y = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = "grey60")
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Fetal kidney reference marker expression by denovo_cell_type"
  )

pdf('fetal.kidney.markers.pdf',height = 7,width=10)
print(p)
dev.off()

#v2

# cleanup
fetal_ <- subset(
  fetal_,
  cells = colnames(fetal_)[!fetal_$denovo_cell_type %in% c(
    "CNT/PC - proximal UB",
    "Pelvic epithelium - distal UB"
  )]
)
fetal_$denovo_cell_type[fetal_$denovo_cell_type %in% "Proximal UB"] <- "UB"

DefaultAssay(fetal_) <- "RNA"
Idents(fetal_) <- "denovo_cell_type"

# ordered cell types within nephron and stroma
nephron_order <- c(
  "Cap mesenchyme",
  "Proliferating cap mesenchyme",
  "Proximal renal vesicle",
  "Distal renal vesicle",
  "Proximal S shaped body",
  "Medial S shaped body",
  "Distal S shaped body",
  "Podocyte",
  "Proximal tubule",
  "Loop of Henle",
  "UB"
)

stroma_order <- c(
  "Stroma progenitor",
  "Proliferating stroma progenitor",
  "Fibroblast",
  "Proliferating fibroblast",
  "Mesangial",
  "Myofibroblast",
  "Smooth muscle"
)

nephron_order <- nephron_order[nephron_order %in% unique(fetal_$denovo_cell_type)]
stroma_order  <- stroma_order[stroma_order %in% unique(fetal_$denovo_cell_type)]

celltype_order <- c(nephron_order, stroma_order)

# markers in blocks
marker_blocks <- list(
  "Cap mesenchyme / blastema-like" = c("SIX2", "CITED1", "EYA1", "WT1", "SALL1", "PAX2"),
  "Early nephron epithelium"       = c("PAX8", "EPCAM", "CDH6"),
  "Tubular"                        = c("LRP2", "SLC12A1"),
  "Podocyte-like"                  = c("NPHS1", "PODXL"),
  "Stromal / smooth muscle"        = c("TAGLN", "MYL9", "PDGFRA"),
  "Myogenic"                       = c("MYOD1", "MYOG", "MYH3")
)

# retain only genes present
marker_blocks <- lapply(marker_blocks, function(x) x[x %in% rownames(fetal_)])
marker_blocks <- marker_blocks[sapply(marker_blocks, length) > 0]

# build Seurat dotplot object
dp <- DotPlot(
  object = fetal_,
  features = marker_blocks,
  group.by = "denovo_cell_type",
  cols = c("lightgrey", "black")
)

plot_df <- dp$data

# order y-axis cell types
plot_df$id <- factor(plot_df$id, levels = rev(celltype_order))

# y-axis grouping
plot_df$component_group <- ifelse(
  as.character(plot_df$id) %in% nephron_order,
  "nephron",
  "stroma"
)
plot_df$component_group <- factor(plot_df$component_group, levels = c("nephron", "stroma"))

# preserve gene order exactly as provided
gene_order <- unlist(marker_blocks, use.names = FALSE)
gene_order <- gene_order[gene_order %in% unique(plot_df$features.plot)]
plot_df$features.plot <- factor(plot_df$features.plot, levels = gene_order)

# optional vertical separator lines between gene blocks
block_lengths <- sapply(marker_blocks, length)
x_breaks <- cumsum(block_lengths)
x_breaks <- x_breaks[-length(x_breaks)] + 0.5

p <- ggplot(plot_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_vline(xintercept = x_breaks, colour = "grey80", linewidth = 0.3) +
  scale_color_gradient(low = "lightgrey", high = "black") +
  facet_grid(
    rows = vars(component_group),
    cols = vars(feature.groups),
    scales = "free",
    space = "free"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = "grey60")
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Scaled expression",
    size = "% cells",
    title = "Fetal kidney marker expression by subtype, grouped into nephron and stroma"
  )

pdf('fetal.kidney.markers.v2.pdf',height = 7,width=10)
print(p)
dev.off()

fetal_ = subset(fetal_, cells = colnames(fetal_)[!fetal_$denovo_cell_type %in% c('CNT/PC - proximal UB','Pelvic epithelium - distal UB')])
fetal_$denovo_cell_type[fetal_$denovo_cell_type %in% "Proximal UB"] = 'UB'

DimPlot(fetal_, group.by = c('denovo_cell_type','component'), label = T, label.box = T, label.size = 2.5, ncol = 1, raster = F, repel = T, shuffle = T) & NoLegend()

nephron <- subset(fetal_,subset=component%in%c('nephron'))
stroma <- subset(fetal_,subset=component%in%c('stroma'))
table(nephron@meta.data$denovo_cell_type)
table(stroma@meta.data$denovo_cell_type)

# Plotting ####

## plotting nephron component ####

# creating a data frame for ggolot
nephron_cells = which(fetal_$component %in% 'nephron')
dt_ = fetal_@meta.data[nephron_cells,]
dt_$x = fetal_[['umap']]@cell.embeddings[nephron_cells,'umap_1']
dt_$y = fetal_[['umap']]@cell.embeddings[nephron_cells,'umap_2']

# finding cluster centers
# labels_ = t(sapply(X = unique(dt_$denovo_cell_type), FUN = function(c_){ return(c(median(dt_$x[dt_$denovo_cell_type %in% c_]),
#                                                                                   median(dt_$y[dt_$denovo_cell_type %in% c_]))) } ))
# colnames(labels_) = c('x','y')
# labels_ = data.frame(labels_, label = rownames(labels_))
label_names <- unique(as.character(dt_$denovo_cell_type))

labels_ <- t(sapply(label_names, function(c_) {
  c(
    x = median(dt_$x[dt_$denovo_cell_type == c_], na.rm = TRUE),
    y = median(dt_$y[dt_$denovo_cell_type == c_], na.rm = TRUE)
  )
}))

labels_ <- data.frame(
  x = labels_[, "x"],
  y = labels_[, "y"],
  label = label_names,
  row.names = NULL
)
# setting cell type sizes to make them distinguishable!

dt_$size = 1
dt_$size[dt_$denovo_cell_type %in% 'Proximal S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Podocyte'] = 4.5

dt_$size[dt_$denovo_cell_type %in% 'Medial S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Proximal tubule'] = 4.5

dt_$size[dt_$denovo_cell_type %in% 'Distal S shaped body'] = 1.5
dt_$size[dt_$denovo_cell_type %in% 'Loop of Henle'] = 4.5


# plotting
ggplot()+
  theme_void()+
  geom_point(data = dt_, aes(x = x, y = y, fill = denovo_cell_type, color = denovo_cell_type), size = dt_$size, show.legend = F, stroke = .2, shape = 21)+
  scale_fill_manual(values = c(
    "Cap mesenchyme"                   = "#deebf7",
    "Proliferating cap mesenchyme"    = "#c6dbef",
    
    "Proximal renal vesicle"          = "#9ecae1",
    "Distal renal vesicle"            = "#6baed6",
    
    "Proximal S shaped body"          = "#4292c6",
    "Medial S shaped body"            = "#3182bd",
    "Distal S shaped body"            = "#2171b5",
    
    "Podocyte"                        = "#084594",
    "Proximal tubule"                 = "#08519c",
    "Loop of Henle"                   = "#08306b",
    
    "UB"                              = "#41b6c4",
    "Stroma progenitor"                 = "#e5f5e0",
    "Proliferating stroma progenitor"  = "#a1d99b",
    "Fibroblast"                        = "#74c476",
    "Proliferating fibroblast"          = "#41ab5d",
    "Mesangial"                         = "#238b45",
    "Myofibroblast"                     = "#006d2c",
    "Smooth muscle"                     = "#00441b"
      
  ))+
  scale_color_manual(values = c("Cap mesenchyme"                   = "#deebf7",
                                "Proliferating cap mesenchyme"    = "#c6dbef",
                                
                                "Proximal renal vesicle"          = "#9ecae1",
                                "Distal renal vesicle"            = "#6baed6",
                                
                                "Proximal S shaped body"          = "#4292c6",
                                "Medial S shaped body"            = "#3182bd",
                                "Distal S shaped body"            = "#2171b5",
                                
                                "Podocyte"                        = "#084594",
                                "Proximal tubule"                 = "#08519c",
                                "Loop of Henle"                   = "#08306b",
                                
                                "UB"                              = "#41b6c4",
                                "Stroma progenitor"                 = "#e5f5e0",
                                "Proliferating stroma progenitor"  = "#a1d99b",
                                "Fibroblast"                        = "#74c476",
                                "Proliferating fibroblast"          = "#41ab5d",
                                "Mesangial"                         = "#238b45",
                                "Myofibroblast"                     = "#006d2c",
                                "Smooth muscle"                     = "#00441b"))+
#geom_text_repel(data = labels_, aes(x = x, y = y, label = label), size = 5, color = "black", fontface = 'bold', nudge_y = .6, seed = 3)
  geom_text_repel(
    data = labels_,
    aes(x = x, y = y, label = label),
    size = 5,
    color = "black",
    fontface = "bold",
    nudge_y = 0.6,
    seed = 3
  )
ggsave(filename = 'fetal_nephron.pdf', device = 'pdf', width = 8, height = 8)

## plotting stroma component ####

## plotting stroma component ####

stroma_cells <- which(fetal_$component %in% "stroma")

dt_stroma <- fetal_@meta.data[stroma_cells, ]
dt_stroma$x <- fetal_[["umap"]]@cell.embeddings[stroma_cells, "umap_1"]
dt_stroma$y <- fetal_[["umap"]]@cell.embeddings[stroma_cells, "umap_2"]

label_names_stroma <- unique(as.character(dt_stroma$denovo_cell_type))

labels_stroma <- t(sapply(label_names_stroma, function(c_) {
  c(
    x = median(dt_stroma$x[dt_stroma$denovo_cell_type == c_], na.rm = TRUE),
    y = median(dt_stroma$y[dt_stroma$denovo_cell_type == c_], na.rm = TRUE)
  )
}))

labels_stroma <- data.frame(
  x = labels_stroma[, "x"],
  y = labels_stroma[, "y"],
  label = label_names_stroma,
  row.names = NULL
)

dt_stroma$size <- 1
dt_stroma$size[dt_stroma$denovo_cell_type %in% "Mesangial"] <- 2.5
dt_stroma$size[dt_stroma$denovo_cell_type %in% "Myofibroblast"] <- 3
dt_stroma$size[dt_stroma$denovo_cell_type %in% "Smooth muscle"] <- 3.5

stroma_cols <- c(
  "Stroma progenitor"                = "#f7fcf5",
  "Proliferating stroma progenitor" = "#d9f0d3",
  "Fibroblast"                       = "#a1d99b",
  "Proliferating fibroblast"         = "#74c476",
  "Mesangial"                        = "#31a354",
  "Myofibroblast"                    = "#006d2c",
  "Smooth muscle"                    = "#00441b"
)

p_stroma <- ggplot() +
  theme_void() +
  geom_point(
    data = dt_stroma,
    aes(x = x, y = y, fill = denovo_cell_type, color = denovo_cell_type),
    size = dt_stroma$size,
    show.legend = FALSE,
    alpha = 0.9,
    stroke = 0.35,
    shape = 21
  ) +
  scale_fill_manual(values = stroma_cols) +
  scale_color_manual(values = stroma_cols) +
  geom_text_repel(
    data = labels_stroma,
    aes(x = x, y = y, label = label),
    size = 5,
    color = "black",
    fontface = "bold",
    nudge_y = 0.6,
    seed = 3
  )

ggsave(
  filename = "fetal_stroma.pdf",
  plot = p_stroma,
  device = "pdf",
  width = 8,
  height = 8
)