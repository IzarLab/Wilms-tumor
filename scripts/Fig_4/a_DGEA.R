# This script finds markers of every malignant cell type of the cancer data (each cell type against all others) and writes the DGEA table
# Run from the root of the repository: Rscript scripts/Fig_4/a_DGEA.R

library(Seurat)
options(Seurat.object.assay.version = "v3")

if(!requireNamespace('presto', quietly = T)){ stop('presto is required: FindAllMarkers uses its fast Wilcoxon test when installed') }

outdir_ = 'out/'
if(!dir.exists(outdir_)){ dir.create(outdir_, recursive = T) }     # create the output directory only if it does not already exist

# STEP 1: Reading in cancer integrated data ####

load('data/compartments/integrated_compartment_cancer_processed.RData')
DefaultAssay(s_objs) = 'RNA'
s_objs[['integrated']] = NULL; gc()     # only the RNA assay is used: frees RAM

# STEP 2: DGEA between all cell types (presto Wilcoxon test) ####

Idents(s_objs) = droplevels(as.factor(s_objs$cell_type))
markers_ = FindAllMarkers(s_objs, only.pos = T, assay = 'RNA', min.pct = 0.07, logfc.threshold = 0.40)
rm(s_objs); gc()     # the Seurat object is not needed anymore

write.csv(markers_, file.path(outdir_, 'DGEA_cancer_cell_types.csv'), row.names = F)
rm(markers_); gc()
