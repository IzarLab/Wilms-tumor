# Wilms Tumor Developmental Plasticity Analysis

## Overview

This repository contains the complete computational workflow used for the analysis of the Wilms tumor developmental plasticity single-nucleus RNA-sequencing study.

The repository includes all preprocessing, dimensional reduction, developmental state annotation, trajectory inference, copy-number analysis, developmental mapping, tumor microenvironment characterization, and figure-generation scripts used throughout the manuscript.

The analyses were performed primarily using:

- R (≥4.3)
- Seurat v5
- Monocle3
- ComplexHeatmap
- circlize
- ggplot2
- patchwork
- pheatmap
- Harmony
- Matrix
- dplyr
- tidyr

---

# Repository structure

```
wilms-tumor-plasticity/

├── data/
│
├── metadata/
│
├── scripts/
│   ├── figure_01/
│   ├── figure_02/
│   ├── figure_03/
│   ├── figure_04/
│   ├── figure_05/
│   └── supplementary/
│
├── figures/
│
├── results/
│
└── README.md
```

---

# Input datasets

The analyses use the following processed Seurat objects.

| Object | Description |
|---------|-------------|
| TME.RData | Integrated tumor microenvironment object |
| integrated_compartment_cancer_processed.RData | Integrated malignant compartment |
| integrated_compartment_lymphoid_processed.RData | Lymphoid compartment |
| integrated_compartment_myeloid_processed.RData | Myeloid compartment |
| integrated_compartment_vasculature_processed.RData | Vasculature compartment |
| Fetal_kidney_processed.RData | Human fetal kidney reference |

Additional metadata tables include

- cohort.tsv
- FGA_TMB.tsv
- clinical metadata
- Numbat CNV outputs
- wmNMF outputs (optional)

---

# Figure 1

## Clinical overview and global atlas

### 01_plot_cohort_circos_summary.R

Generates the manuscript cohort summary circos plot integrating

- clinical characteristics
- tumor cellularity
- FGA
- TMB
- treatment
- histology
- response
- recurrence
- survival
- tumor laterality
- tumor size
- stage
- neoadjuvant regimen

Outputs

- Cohort circos plot
- Legend PDF

---

### 02_plot_global_umaps.R

Generates UMAP visualizations of the complete integrated dataset.

Panels include

- major compartments
- treatment status
- cell-cycle phase

Outputs

- labeled and unlabeled UMAPs

---

### 03_plot_malignant_umaps_composition_markers.R

Analyzes the malignant compartment.

Generates

- malignant UMAPs
- treatment UMAPs
- stacked composition plots
- dot plots
- differential abundance summaries

---

# Figure 2

## Developmental programs

### 01_plot_malignant_subcompartment_analysis.R

Produces

- malignant UMAPs
- condition comparisons
- stacked bar plots
- abundance statistics

---

### 02_plot_fetal_kidney_reference.R

Maps Wilms tumor populations onto the human fetal kidney atlas.

Generates

- fetal nephron marker dot plots
- stromal marker dot plots
- developmental reference figures

---

### 03_plot_nmf_metaprogram_heatmap.R

Constructs the cleaned wmNMF metaprogram heatmap.

Inputs

- Jaccard similarity matrix
- annotated metaprograms

Outputs

- manuscript-ready metaprogram heatmap

---

### 04_plot_top_marker_dotplots.R

Computes

- marker genes
- top marker dot plots

for

- lymphoid
- myeloid
- vasculature

---

# Figure 3

## Developmental plasticity

### 01_monocle_untreated_plasticity.R

Constructs developmental trajectories for untreated tumors.

Includes

- Monocle3 trajectory inference
- pseudotime estimation
- developmental transition graph
- transition flux analysis

---

### 02_cnv_plasticity_summary.R

Summarizes inferred copy-number profiles.

Produces

- chromosome-arm CNV heatmaps
- CNV burden comparisons
- compartment-level CNV statistics

---

### 03_patient_candle_summary.R

Generates patient-level developmental plasticity summaries using the CANDLE framework.

Outputs

- per-patient switching metrics
- within-state plasticity
- between-state selection
- patient summary panels

---

# Figure 4

## Treatment-associated developmental remodeling

### 01_monocle_joint_trajectory.R

Constructs a joint developmental trajectory containing

- untreated tumors
- treated tumors

Generates

- shared trajectory graph
- pseudotime
- transition networks

---

### 02_monocle_treatment_regimens.R

Builds treatment-specific developmental trajectories.

Regimens

- EE-4A
- DD-4A
- Other

Outputs

- regimen-specific transition graphs
- developmental flux
- trajectory comparisons

---

### 03_treatment_response_analysis.R

Analyzes treatment-associated developmental remodeling.

Includes

- responder/non-responder comparisons
- developmental category analysis
- clinical metadata integration

---

# Figure 5

## Tumor microenvironment

### 01_plot_tme_composition.R

Analyzes immune and stromal composition.

Includes

- lymphoid
- myeloid
- vasculature

Outputs

- stacked composition plots
- clinical annotation panels
- treatment comparisons

---

# Software requirements

Required R packages include

```
Seurat
SeuratObject
monocle3
ggplot2
patchwork
ComplexHeatmap
circlize
pheatmap
Harmony
Matrix
dplyr
tidyr
tibble
readr
scales
igraph
ggraph
cowplot
viridis
rstatix
```

---

# Running the analysis

The scripts are designed to be executed independently after loading the required processed Seurat objects.

Typical workflow:

```
Figure 1
↓

Figure 2
↓

Figure 3
↓

Figure 4
↓

Figure 5
```

Intermediate processed objects generated during one analysis can be reused by subsequent scripts.

---

# Notes

- All scripts were used to generate the figures included in the manuscript.

- File paths should be updated to match the local project directory.

- Output directories are automatically created when missing.

- Random seeds are fixed where stochastic algorithms are used to ensure reproducibility.

---

# Citation

To be updated.
