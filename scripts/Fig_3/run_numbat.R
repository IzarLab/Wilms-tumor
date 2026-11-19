#!/usr/bin/env Rscript

# Run haplotype-aware CNV inference for one Wilms tumor sample.
#
# Please change the path and provide associated input files as required. The provided files are just for reference.
#
# Typical Usage from the repository root
#   Rscript scripts/Fig_3/run_numbat.R \
#     inputs/W001_counts.rds inputs/W001_allele_counts.tsv.gz \
#     inputs/lambdas_ref.rds inputs/numbat_options.rds numbat_output/W001
#
# Dependencies
#   R packages: numbat, data.table, Matrix (for sparse count matrices).
#   Use the Numbat release recorded for the analysis you are reproducing.
#
# Inputs, in positional order
#   1. COUNTS.rds: a genes-by-nuclei raw RNA UMI count matrix, preferably a
#      sparse dgCMatrix, with gene IDs as row names and unique barcodes as
#      column names. Use the RNA counts layer before normalization.
#   2. ALLELES.tsv.gz: Numbat's phased allele-count table, with at least the
#      documented fields cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT
#      and gene. Use the same genome build for the alignment, genetic map,
#      SNP VCF, phasing panel and Numbat configuration.
#   3. LAMBDAS_REF.rds: the normalized expression reference supplied to Numbat.
#      Use the reference from the analysis being reproduced. For a custom
#      reference, aggregate raw counts by annotated reference group with
#      numbat::aggregate_counts(ref_counts, ref_annotation), where
#      ref_annotation contains columns cell and group. Save the resulting
#      reference matrix with saveRDS(). Alternatively, explicitly save a
#      package-provided reference if that is the intended input.
#   4. OPTIONS.rds: a named R list of run_numbat() arguments. Explicitly
#      specify genome, max_iter, min_LLR, gamma and t, plus any overrides.
#   5. OUT_DIR: a new or empty sample directory for Numbat outputs.
#
# Upstream allele preparation
#   Use the pileup_and_phase.R script from the installed Numbat release. The
#   standard workflow uses cellsnp-lite for SNP pileup and Eagle2 for phasing.
#   Run Rscript /path/to/pileup_and_phase.R --help for release-specific flags.
#   Typical arguments include --label, --samples, --bams, --barcodes, --gmap,
#   --eagle, --snpvcf, --paneldir, --outdir and --ncores. Retain the full
#   preprocessing command and record tool versions and reference identifiers.
#   See https://kharchenkolab.github.io/numbat/articles/numbat.html
#
# Example options for a new run using Numbat v1.5.2 defaults
#   saveRDS(list(genome = "hg38", max_iter = 2L, min_LLR = 5,
#                gamma = 20, t = 1e-5), "inputs/numbat_options.rds")
#   This is an example configuration. Historical Wilms preprocessing,
#   expression reference, package version and run settings require the
#   original run records. The first-iteration Fig. 3a input filename does
#   not identify the configured maximum number of iterations.
#
# Nucleus selection and outputs
#   The inference uses all barcodes shared by the raw counts and phased
#   allele table. Malignancy annotation is applied downstream. The script
#   records the run options and session alongside Numbat's results. Fig. 3a
#   consumes joint_post_1.tsv through summarize_cnv.R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript run_numbat.R COUNTS.rds ALLELES.tsv.gz LAMBDAS_REF.rds OPTIONS.rds OUT_DIR")
}

suppressPackageStartupMessages({
  library(data.table)
  library(numbat)
})

counts <- readRDS(args[1])
alleles <- fread(args[2], colClasses = c(cell = "character"))
reference <- readRDS(args[3])
run_options <- readRDS(args[4])
output <- args[5]

# Keep data inputs and output location under the control of the CLI.
required <- c("genome", "max_iter", "min_LLR", "gamma", "t")
reserved <- c("count_mat", "df_allele", "lambdas_ref", "out_dir")
option_names <- names(run_options)
if (!is.list(run_options) || anyNA(option_names) ||
    any(!nzchar(option_names)) || anyDuplicated(option_names) ||
    !all(required %in% option_names) || any(option_names %in% reserved)) {
  stop("OPTIONS.rds must be a uniquely named list containing genome, max_iter, min_LLR, gamma and t.")
}

unsupported <- setdiff(option_names, names(formals(numbat::run_numbat)))
if (length(unsupported)) {
  stop("Options unavailable in this Numbat release: ", paste(unsupported, collapse = ", "))
}

# Match exactly. Apply any sequencing-pipeline barcode normalization to the
# inputs before running this script, so different nuclei cannot be merged.
if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
    anyNA(colnames(counts)) || anyDuplicated(colnames(counts))) {
  stop("The count matrix must have gene row names and unique, nonmissing barcode column names.")
}
if (!"cell" %in% names(alleles)) {
  stop("The phased allele table must contain a cell column.")
}
cells <- intersect(colnames(counts), unique(alleles$cell))
if (!length(cells)) {
  stop("No shared barcodes between the count matrix and phased allele table.")
}
message("Using ", length(cells), " of ", ncol(counts), " expression barcodes.")

# Use a clean directory to avoid mixing posterior tables from separate runs.
if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
  stop("Choose a new or empty output directory: ", output)
}
if (!dir.exists(output) && !dir.create(output, recursive = TRUE)) {
  stop("Could not create output directory: ", output)
}
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
saveRDS(run_options, file.path(output, "numbat_options.rds"))

call_args <- c(list(
  count_mat = counts[, cells, drop = FALSE],
  df_allele = alleles[cell %in% cells],
  lambdas_ref = reference,
  out_dir = output
), run_options)
do.call(numbat::run_numbat, call_args)
