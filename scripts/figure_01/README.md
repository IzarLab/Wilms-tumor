# Figure 1 scripts

This folder contains GitHub-ready scripts used to generate Figure 1 panels for the Wilms tumor plasticity manuscript.

## Folder placement

Place this folder inside the main repository as:

```text
scripts/figure_01/
```

Recommended repository structure:

```text
wilms-tumor-plasticity/
├── scripts/
│   └── figure_01/
├── data/
│   ├── metadata/
│   └── processed/
├── results/
│   └── figure_01/
└── figures/
    └── figure_01/
```

## Scripts

### `01_plot_cohort_circos_summary.R`

Generates the cohort-level circos summary plot integrating malignant subcompartment composition, FGA, TMB, demographics, clinical metadata, pathology, treatment, neoadjuvant regimen, and response.

Required inputs:

```text
data/metadata/cohort.tsv
data/metadata/FGA_TMB.tsv
data/processed/TME.RData
```

Expected object inside `TME.RData`:

```r
tme_
```

Outputs:

```text
figures/figure_01/cohort_summary_circos.pdf
figures/figure_01/cohort_summary_circos_legend_only.pdf
```

Run:

```bash
Rscript scripts/figure_01/01_plot_cohort_circos_summary.R
```

### `02_plot_global_umaps.R`

Generates global TME UMAPs colored by major compartment, treatment condition, and cell-cycle phase.

Required input:

```text
data/processed/TME.RData
```

Expected object inside `TME.RData`:

```r
tme_
```

Outputs:

```text
figures/figure_01/global_umap_compartment_no_label.pdf
figures/figure_01/global_umap_compartment_label.pdf
figures/figure_01/global_umap_condition_no_label.pdf
figures/figure_01/global_umap_condition_label.pdf
figures/figure_01/global_umap_phase_no_label.pdf
figures/figure_01/global_umap_phase_label.pdf
```

Run:

```bash
Rscript scripts/figure_01/02_plot_global_umaps.R
```

### `03_plot_malignant_umaps_composition_markers.R`

Generates malignant-cell UMAPs, malignant subcompartment composition plots, top-marker dot plots, and per-sample abundance statistics.

Required input:

```text
data/processed/integrated_compartment_cancer_processed.RData
```

Expected object inside the RData file:

```r
s_objs
```

Outputs:

```text
figures/figure_01/malignant_umap_subcompartment_no_label.pdf
figures/figure_01/malignant_umap_subcompartment_label.pdf
figures/figure_01/malignant_umap_condition_no_label.pdf
figures/figure_01/malignant_umap_condition_label.pdf
figures/figure_01/malignant_stacked_bar_condition_by_subcompartment.pdf
figures/figure_01/malignant_stacked_bar_condition_by_subcompartment_labeled.pdf
figures/figure_01/malignant_dotplot_top10_markers_per_subcompartment.pdf
figures/figure_01/malignant_boxplot_subcompartment_fraction_by_condition.pdf
results/figure_01/all_markers_subcompartment.csv
results/figure_01/top10_markers_per_subcompartment.csv
results/figure_01/subcompartment_fraction_by_condition_stats.csv
```

Run:

```bash
Rscript scripts/figure_01/03_plot_malignant_umaps_composition_markers.R
```

## Notes

All scripts use relative paths from the repository root. Run them from the root directory of the GitHub repository, not from inside `scripts/figure_01/`.
