#!/usr/bin/env Rscript
# R Detectors: scDblFinder (native Call with dbr.sd=1), cxds, bcds, hybrid,
# and DoubletFinder pANN scores. DoubletFinder classification is discarded.

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

parse_args <- function(args) {
  out <- list(mtx_dir = NULL, outdir = NULL, doubletfinder = FALSE)
  i <- 1
  while (i <= length(args)) {
    if (args[[i]] == "--mtx-dir") {
      out$mtx_dir <- args[[i + 1]]
      i <- i + 2
    } else if (args[[i]] == "--outdir") {
      out$outdir <- args[[i + 1]]
      i <- i + 2
    } else if (args[[i]] == "--doubletfinder") {
      out$doubletfinder <- TRUE
      i <- i + 1
    } else {
      i <- i + 1
    }
  }
  if (is.null(out$mtx_dir) || is.null(out$outdir)) {
    stop("usage: detectors_r.R --mtx-dir DIR --outdir DIR [--doubletfinder]")
  }
  out
}

write_skip <- function(outdir, name, reason) {
  writeLines(reason, file.path(outdir, paste0(name, ".skip.txt")))
}

write_tsv <- function(outdir, name, barcodes, score, native_call = NULL) {
  df <- data.frame(barcode = as.character(barcodes), score = as.numeric(score), stringsAsFactors = FALSE)
  if (!is.null(native_call)) {
    df$native_call <- as.character(native_call)
  }
  utils::write.table(df, file.path(outdir, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

load_counts <- function(mtx_dir) {
  mat <- Matrix::readMM(file.path(mtx_dir, "matrix.mtx"))
  barcodes <- readLines(file.path(mtx_dir, "barcodes.tsv"))
  genes <- readLines(file.path(mtx_dir, "genes.tsv"))
  mat <- as(mat, "dgCMatrix")
  if (length(genes) != nrow(mat)) stop("genes.tsv length != nrow(matrix)")
  if (length(barcodes) != ncol(mat)) stop("barcodes.tsv length != ncol(matrix)")
  rownames(mat) <- make.unique(genes)
  colnames(mat) <- make.unique(barcodes)
  mat
}

opt <- parse_args(args)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(0)

mat <- tryCatch(load_counts(opt$mtx_dir), error = function(e) {
  write_skip(opt$outdir, "scdblfinder", paste0("error:", conditionMessage(e)))
  write_skip(opt$outdir, "cxds", paste0("error:", conditionMessage(e)))
  write_skip(opt$outdir, "bcds", paste0("error:", conditionMessage(e)))
  write_skip(opt$outdir, "hybrid", paste0("error:", conditionMessage(e)))
  if (opt$doubletfinder) write_skip(opt$outdir, "doubletfinder", paste0("error:", conditionMessage(e)))
  quit(save = "no", status = 0)
})

barcodes <- colnames(mat)

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  write_skip(opt$outdir, "scdblfinder", "missing_package:SingleCellExperiment")
  write_skip(opt$outdir, "cxds", "missing_package:SingleCellExperiment")
  write_skip(opt$outdir, "bcds", "missing_package:SingleCellExperiment")
  write_skip(opt$outdir, "hybrid", "missing_package:SingleCellExperiment")
} else {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))

  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    write_skip(opt$outdir, "scdblfinder", "missing_package:scDblFinder")
  } else {
    tryCatch({
      # dbr.sd=1: threshold from artificial-doublet misclassification, not Expected Doublet Rate.
      dbl <- scDblFinder::scDblFinder(sce, clusters = NULL, dbr.sd = 1, verbose = TRUE)
      write_tsv(
        opt$outdir, "scdblfinder", colnames(dbl),
        dbl$scDblFinder.score,
        as.character(dbl$scDblFinder.class)
      )
    }, error = function(e) write_skip(opt$outdir, "scdblfinder", paste0("error:", conditionMessage(e))))
  }

  if (!requireNamespace("scds", quietly = TRUE)) {
    write_skip(opt$outdir, "cxds", "missing_package:scds")
    write_skip(opt$outdir, "bcds", "missing_package:scds")
    write_skip(opt$outdir, "hybrid", "missing_package:scds")
  } else {
    tryCatch({
      sce_c <- scds::cxds(sce)
      write_tsv(opt$outdir, "cxds", colnames(sce_c), sce_c$cxds_score)
    }, error = function(e) write_skip(opt$outdir, "cxds", paste0("error:", conditionMessage(e))))
    tryCatch({
      sce_b <- scds::bcds(sce, verb = FALSE)
      write_tsv(opt$outdir, "bcds", colnames(sce_b), sce_b$bcds_score)
    }, error = function(e) write_skip(opt$outdir, "bcds", paste0("error:", conditionMessage(e))))
    tryCatch({
      sce_h <- scds::cxds_bcds_hybrid(sce)
      write_tsv(opt$outdir, "hybrid", colnames(sce_h), sce_h$hybrid_score)
    }, error = function(e) write_skip(opt$outdir, "hybrid", paste0("error:", conditionMessage(e))))
  }
}

if (isTRUE(opt$doubletfinder)) {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("DoubletFinder", quietly = TRUE)) {
    missing <- c()
    if (!requireNamespace("Seurat", quietly = TRUE)) missing <- c(missing, "Seurat")
    if (!requireNamespace("DoubletFinder", quietly = TRUE)) missing <- c(missing, "DoubletFinder")
    write_skip(opt$outdir, "doubletfinder", paste0("missing_package:", paste(missing, collapse = ",")))
  } else {
    tryCatch({
      seu <- Seurat::CreateSeuratObject(counts = mat)
      seu <- Seurat::NormalizeData(seu, verbose = FALSE)
      seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      npcs <- min(10L, ncol(seu) - 1L)
      seu <- Seurat::RunPCA(seu, npcs = npcs, verbose = FALSE)
      sweep.res <- tryCatch(
        DoubletFinder::paramSweep(seu, PCs = 1:npcs, sct = FALSE),
        error = function(e) DoubletFinder::paramSweep_v3(seu, PCs = 1:npcs, sct = FALSE)
      )
      sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
      bcmvn <- DoubletFinder::find.pK(sweep.stats)
      pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      # nExp=1 satisfies the API only. The DF class column is discarded; Call is MAD on pANN.
      seu <- DoubletFinder::doubletFinder(seu, PCs = 1:npcs, pN = 0.25, pK = pK, nExp = 1L)
      md <- seu@meta.data
      pann_col <- grep("^pANN_", colnames(md), value = TRUE)
      if (length(pann_col) == 0) stop("no pANN column")
      write_tsv(opt$outdir, "doubletfinder", colnames(seu), md[[pann_col[[1]]]])
    }, error = function(e) write_skip(opt$outdir, "doubletfinder", paste0("error:", conditionMessage(e))))
  }
}
