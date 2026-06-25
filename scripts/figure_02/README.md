# Figure 2 scripts

Scripts for malignant compartment annotation, developmental-state marker visualization, fetal kidney reference marker comparison, and NMF metaprogram heatmap generation.

Recommended run order:
1. `01_plot_malignant_umaps_composition_markers.R`
2. `02_plot_nmf_program_jaccard_heatmap.R`
3. `03_plot_fetal_kidney_reference_markers.R`

Expected major inputs include the processed malignant Seurat object, NMF Jaccard matrix, manually curated NMF program annotation table, and fetal kidney reference Seurat object.
