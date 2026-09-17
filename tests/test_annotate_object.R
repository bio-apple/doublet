#!/usr/bin/env Rscript
# Unit checks for Seurat / SingleCellExperiment column attach (no Detector run).

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
} else {
  getwd()
}
source(file.path(root, "R", "annotate_doublets.R"))

cells <- data.frame(
  barcode = c("c1", "c2", "c3"),
  scdblfinder_score = c(0.1, 0.9, 0.2),
  scdblfinder_call = c("singlet", "doublet", "singlet"),
  scrublet_score = c(0.4, 0.5, 0.6),
  scrublet_call = c("singlet", "singlet", "doublet"),
  stringsAsFactors = FALSE
)
sample <- data.frame(
  detector = c("scdblfinder", "scrublet", "solo"),
  status = c("ran", "ran", "skipped"),
  stringsAsFactors = FALSE
)

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  stop("SingleCellExperiment required for this test")
}
mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
sce <- attach_doublet_columns(sce, cells, sample, primary = "scdblfinder")
stopifnot(isTRUE(all.equal(as.numeric(sce$doublet_score), c(0.1, 0.9, 0.2))))
stopifnot(identical(as.logical(sce$predicted_doublet), c(FALSE, TRUE, FALSE)))
stopifnot(identical(as.logical(sce$is_doublet), c(FALSE, TRUE, FALSE)))
stopifnot(identical(S4Vectors::metadata(sce)$doublet_primary, "scdblfinder"))
stopifnot(identical(counts_from_object(sce), SingleCellExperiment::counts(sce)))

if (requireNamespace("Seurat", quietly = TRUE)) {
  seu <- Seurat::CreateSeuratObject(counts = mat)
  seu <- attach_doublet_columns(seu, cells, sample, primary = "scdblfinder")
  stopifnot(isTRUE(all.equal(as.numeric(seu$doublet_score), c(0.1, 0.9, 0.2))))
  stopifnot(identical(as.logical(seu$predicted_doublet), c(FALSE, TRUE, FALSE)))
  stopifnot(identical(as.logical(seu$is_doublet), c(FALSE, TRUE, FALSE)))
}

cat("ok\n")
