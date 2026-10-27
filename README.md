# Wilms Tumor Developmental Plasticity Analysis

R scripts that generate the figures of the Wilms tumor developmental plasticity single-nucleus RNA-sequencing study. Each script
produces one figure panel (or a few related panels) from processed Seurat objects and tables in `data/`.

## Repository structure

```
wilms-tumor-plasticity/
├── data/
│   ├── compartments/     # processed Seurat objects (.RData)
│   └── Misc/             # tables, per-sample NMF results and Milo objects
├── scripts/
│   ├── Fig_1/ ... Fig_5/  # main figures
│   └── ExtFig_1/, ExtFig_3/  # Extended Data figures
├── out/                  # created on first run; all figures are written here
└── README.md
```

Script names follow the panels they generate, e.g. `scripts/Fig_2/b-c.R` makes Fig. 2b and 2c, and `scripts/Fig_5/b_d_f.R`
makes Fig. 5b, 5d and 5f.

## Input data

Seurat objects in `data/compartments/`:

| File | Content |
|------|---------|
| `TME.RData` | Integrated tumor microenvironment (all cells) |
| `integrated_compartment_cancer_processed.RData` | Integrated malignant compartment |
| `integrated_compartment_lymphoid_processed.RData` | Integrated lymphoid compartment |
| `integrated_compartment_myeloid_processed.RData` | Integrated myeloid compartment |
| `integrated_compartment_vasculature_processed.RData` | Integrated vasculature compartment |
| `Fetal_kideny_processed.RData` | Human fetal kidney reference (kidney cell atlas, Stewart et al. 2019) |

Tables and intermediate results in `data/Misc/`:

| File | Content |
|------|---------|
| `cohort.tsv` | Patient clinical metadata, including neoadjuvant regimen |
| `FGA_TMB.tsv` | Fraction of genome altered and tumor mutational burden per patient |
| `ALL_samples_*.tsv` | Copy-number (CNV) tables per sample, chromosome arm and cell |
| `NMF_W*.RData` | Per-sample NMF results |
| `Jaccard_similarity_matrix_programs.tsv` | Jaccard similarity between NMF factors |
| `programs_annotation_processed.tsv` | Manually annotated NMF programs (factor-sample pairs) |
| `celltypes_markers_descartes.tsv`, `celltypes_markers_kidneycellatlas.tsv` | Fetal kidney cell type markers |
| `Milo_DA_k30_d30_prop0.4_*.RData` | Milo neighborhood objects (cancer, lymphoid, myeloid, vasculature) |

## How to run

Every script is run **from the root of the repository** with `Rscript`, for example:

```bash
Rscript scripts/Fig_2/b-c.R
```

- Scripts are independent of each other and can be run in any order.
- Figures are written to `out/`. The folder is created if it does not exist; existing files with the same name are overwritten.
- Random seeds are fixed in scripts with stochastic steps (downsampling, Milo, Monocle3), so results are reproducible.
- The Seurat objects are large. Scripts subset them and remove them from memory as soon as they are no longer needed, but
  scripts that load the cancer or TME object still need a machine with plenty of RAM.

### Scripts

| Script | Panel(s) | What it does | Inputs | Outputs (in `out/`) |
|--------|----------|--------------|--------|---------------------|
| `Fig_1/b.R` | 1b | Circos summary of the cohort: clinical, pathology, genomic (FGA, TMB) and subcompartment features per patient | `TME.RData`, `cohort.tsv`, `FGA_TMB.tsv` | `cohort_summary_circos.pdf`, `cohort_summary_circos_legend_only.pdf` |
| `Fig_1/c-e.R` | 1c-e | UMAPs of the whole TME colored by patient, condition, compartment, subcompartment, cell type and cell-cycle phase | `TME.RData` | `global_umaps.pdf` |
| `Fig_2/a.R` | 2a | UMAP of the malignant compartment colored by subcompartment | cancer object | `WT_component_umap.pdf` |
| `Fig_2/b-c.R` | 2b, 2c | Top 10 protein-coding markers per malignant subcompartment (dot plot) and subcompartment abundance by condition (Wilcoxon, BH) | cancer object | `dotplot_top10_markers_per_subcompartment.pdf`, `boxplot_subcompartment_fraction_by_condition_with_pvalues.pdf`, `all_markers_subcompartment.csv`, `top10_markers_per_subcompartment.csv` |
| `Fig_2/d.R` | 2d | Jaccard similarity heatmap of NMF factors clustered into programs | `Jaccard_similarity_matrix_programs.tsv`, `programs_annotation_processed.tsv` | `SOM_NMF_programs_cleaned.pdf` |
| `Fig_2/e.R` | 2e | Top 5 genes of each NMF program; color = mean expression in expressing cells of the program's cell type, size = NMF weight | cancer object, `NMF_W*.RData`, `programs_annotation_processed.tsv` | `NMF_program_markers_dotplot.pdf` |
| `Fig_2/f.R` | 2f | Maps NMF programs to fetal kidney cell types by overlap of program genes and fetal markers (overlap coefficient, hypergeometric test) | fetal object, `NMF_W*.RData`, `programs_annotation_processed.tsv` | `NMF_programs_fetal_kidney_overlap_dotplot.pdf` |
| `Fig_3/a.R` | 3a | Chromosome-arm CNV state and CNV burden per subcompartment in untreated samples | `ALL_samples_*.tsv` | `Wilms_UNTREATED_compact_PanelE_left_PanelC_right.pdf/.png` |
| `Fig_3/b-c.R` | 3b, 3c | Monocle3 trajectory of untreated malignant cells, rooted at the most potent cell type (CytoTRACE 2), and pseudotime by subcompartment | cancer object | `UMAP_monocle3_trajectory_untreated.pdf`, `boxplot_pseudotime_by_subcompartment.pdf`, `CytoTRACE2_pct_potent_per_celltype.csv` |
| `Fig_4/a.R` | 4a | UMAP of the malignant compartment colored by cell type | cancer object | `malignant_compartment_umap.pdf/.png` |
| `Fig_4/b.R` | 4b | Milo differential abundance of malignant cell neighborhoods, neoadjuvant vs untreated | cancer object, `Milo_DA_k30_d30_prop0.4_cancer.RData` | `Milo_DA_k30_d30_prop0.4_cancer_*_UMAP.pdf`, `Milo_DA_k30_d30_prop0.4_cancer_stats.pdf` |
| `Fig_4/c.R` | 4c | Ternary density plot of treated malignant cells by fetal blastema, epithelium and stroma signature scores | cancer object, fetal marker tables | `cancer_density_plasticity.pdf`, `cancer_density_plasticity_subcompartments.pdf` |
| `Fig_4/d.R` | 4d | Joint Monocle3 trajectory of untreated and treated cells, colored by treatment enrichment | cancer object | `joint_trajectory_condition_enrichment.pdf`, `joint_trajectory_celltype_enrichment.csv`, `CytoTRACE2_pct_potent_per_celltype_joint.csv` |
| `Fig_4/e.R` | 4e | Monocle3 trajectory of treated cells with the share of each regimen (EE-4A, DD-4A, other) per cell type | cancer object, `cohort.tsv` | `treated_trajectory_regimen.pdf`, `treated_trajectory_celltype_regimen.csv`, `CytoTRACE2_pct_potent_per_celltype_treated.csv` |
| `Fig_5/a.R` | 5a | UMAP of the nonmalignant TME colored by subcompartment | `TME.RData` | `nonmalignant_TME.pdf` |
| `Fig_5/b_d_f.R` | 5b, 5d, 5f | UMAPs of the lymphoid, myeloid and vasculature compartments colored by cell type | lymphoid, myeloid and vasculature objects | `lymphoid_compartment_umap.pdf`, `myeloid_compartment_umap.pdf`, `vasculature_compartment_umap.pdf` |
| `Fig_5/c.R` | 5c | Milo differential abundance of lymphoid neighborhoods | lymphoid object, `Milo_DA_k30_d30_prop0.4_lymphoid.RData` | `Milo_DA_k30_d30_prop0.4_lymphoid_*` |
| `Fig_5/e.R` | 5e | Milo differential abundance of myeloid neighborhoods | myeloid object, `Milo_DA_k30_d30_prop0.4_myeloid.RData` | `Milo_DA_k30_d30_prop0.4_myeloid_*` |
| `Fig_5/g.R` | 5g | Milo differential abundance of vascular neighborhoods | vasculature object, `Milo_DA_k30_d30_prop0.4_vasculature.RData` | `Milo_DA_k30_d30_prop0.4_vasculature_*` |
| `ExtFig_1/a.R` | ED 1a | Canonical fetal markers selected from a DGEA between malignant cell types, plotted across fetal kidney cell types | cancer and fetal objects | `fetal_markers_dotplot.pdf` |
| `ExtFig_1/b.R` | ED 1b | Maps malignant cell types to fetal kidney cell types by marker overlap (overlap coefficient, hypergeometric test) | cancer and fetal objects | `malignant_cell_types_fetal_kidney_overlap_dotplot.pdf` |
| `ExtFig_1/f-g.R` | ED 1f, 1g | UMAPs of the fetal kidney nephron and stroma components | fetal object | `fetal_nephron.pdf`, `fetal_stroma.pdf` |
| `ExtFig_3/a.R` | ED 3a | Scaled pseudotime per cell type and condition (joint trajectory of Fig. 4d) | cancer object | `joint_pseudotime_by_celltype_condition.pdf` |
| `ExtFig_3/b.R` | ED 3b | Scaled pseudotime per cell type and regimen (treated trajectory of Fig. 4e) | cancer object, `cohort.tsv` | `treated_pseudotime_by_celltype_regimen.pdf` |

"cancer object" is `integrated_compartment_cancer_processed.RData` and "fetal object" is `Fetal_kideny_processed.RData`.

### Notes

- **Milo scripts** (`Fig_4/b.R`, `Fig_5/c.R`, `Fig_5/e.R`, `Fig_5/g.R`) read the Milo object from `data/Misc/` if it exists.
  Otherwise they build the neighborhoods and save the object there, which takes much longer.
- **`Fig_2/4_NMF_program_heatmap.R`** is an upstream step, not a figure script. It builds the Jaccard similarity matrix and the
  unannotated program table (`programs_annotation.tsv`) from the per-sample NMF and GSEA outputs, which are not part of this
  repository. Its annotated output, `programs_annotation_processed.tsv`, is provided in `data/Misc/`.
- **CytoTRACE 2** is called as `CytoTRACE2::cytotrace2()` and not attached, because attaching it masks `dplyr` verbs.
- Scripts ending in `_old.R` are earlier versions kept for reference; they are not needed to reproduce the figures.

## Software

All scripts were run with **R 4.5.3** (2026-03-11) and **Bioconductor 3.22** on macOS (aarch64-apple-darwin20).

| Package | Version | Source |
|---------|---------|--------|
| AnnotationDbi | 1.72.0 | Bioconductor |
| circlize | 0.4.18 | CRAN |
| ComplexHeatmap | 2.26.1 | Bioconductor |
| CytoTRACE2 | 1.1.0 | GitHub (digitalcytometry/cytotrace2) |
| dplyr | 1.2.1 | CRAN |
| edgeR | 4.8.2 | Bioconductor (used by miloR) |
| FNN | 1.1.4.1 | CRAN |
| ggforce | 0.5.0 | CRAN |
| ggplot2 | 4.0.3 | CRAN |
| ggpubr | 1.0.0 | CRAN |
| ggrepel | 0.9.8 | CRAN |
| ggtern | 4.0.0 | CRAN |
| gridExtra | 2.3.1 | CRAN |
| igraph | 2.3.4 | CRAN |
| Matrix | 1.7.6 | CRAN |
| miloR | 2.6.0 | Bioconductor |
| monocle3 | 1.4.27 | GitHub (cole-trapnell-lab/monocle3) |
| NMF | 0.28 | CRAN |
| org.Hs.eg.db | 3.22.0 | Bioconductor |
| patchwork | 1.3.2 | CRAN |
| pheatmap | 1.0.13 | CRAN |
| presto | 1.1.0 | GitHub (immunogenomics/presto; used by Seurat's `FindAllMarkers`) |
| readr | 2.2.0 | CRAN |
| rstatix | 1.1.0 | CRAN |
| scales | 1.4.0 | CRAN |
| scater | 1.38.1 | Bioconductor |
| scatterpie | 0.2.6 | CRAN |
| Seurat | 5.6.0 | CRAN |
| SeuratObject | 5.4.0 | CRAN |
| SingleCellExperiment | 1.32.0 | Bioconductor |
| statmod | 1.5.2 | CRAN (used by edgeR) |
| tidyr | 1.3.2 | CRAN |
| tidyverse | 2.0.0 | CRAN |
| viridis | 0.6.5 | CRAN |

Base R packages used: `grid`, `tools`.

### Installation

```r
install.packages(c('BiocManager', 'remotes'))
install.packages(c('circlize', 'dplyr', 'FNN', 'ggforce', 'ggplot2', 'ggpubr', 'ggrepel', 'ggtern', 'gridExtra', 'igraph',
                   'Matrix', 'NMF', 'patchwork', 'pheatmap', 'readr', 'rstatix', 'scales', 'scatterpie', 'Seurat',
                   'SeuratObject', 'statmod', 'tidyr', 'tidyverse', 'viridis'))
BiocManager::install(c('AnnotationDbi', 'ComplexHeatmap', 'edgeR', 'miloR', 'org.Hs.eg.db', 'scater', 'SingleCellExperiment'))
remotes::install_github('cole-trapnell-lab/monocle3')
remotes::install_github('immunogenomics/presto')
remotes::install_github('digitalcytometry/cytotrace2', subdir = 'cytotrace2_r')
```

## Citation

If you use this code, please cite:

Chalabi Hajkarim M. et al. Cancer cell plasticity through the lens of Wilms tumor developmental state dynamics. (Manuscript in preparation)
