#!/usr/bin/env Rscript
# Annotate Seurat or SingleCellExperiment with Detector Scores/Calls.
# source() this file, then: x <- annotate_doublets(x, fast = TRUE)

DETECTORS <- c(
  "scdblfinder", "scrublet", "cxds", "bcds", "hybrid",
  "doubletdetection", "doubletfinder", "solo"
)

annotate_script_dir <- function() {
  ofile <- NULL
  n <- sys.nframe()
  if (n >= 1) {
    for (i in n:1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) break
    }
  }
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile)))
  }
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath("scripts")
}

counts_from_object <- function(x) {
  if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("package SingleCellExperiment is required")
    }
    names <- SummarizedExperiment::assayNames(x)
    if (!("counts" %in% names)) {
      stop("SingleCellExperiment has no 'counts' assay")
    }
    return(SingleCellExperiment::counts(x))
  }
  if (inherits(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("package Seurat is required")
    }
    assay <- Seurat::DefaultAssay(x)
    mat <- tryCatch(
      Seurat::GetAssayData(object = x, assay = assay, layer = "counts"),
      error = function(e) Seurat::GetAssayData(object = x, assay = assay, slot = "counts")
    )
    if (is.list(mat)) {
      mat <- mat[[1]]
    }
    return(mat)
  }
  stop("expected a Seurat or SingleCellExperiment object")
}

write_counts_mtx <- function(mat, mtx_dir) {
  dir.create(mtx_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(colnames(mat)) || is.null(rownames(mat))) {
    stop("count matrix needs rownames (genes) and colnames (barcodes)")
  }
  mat <- as(mat, "dgCMatrix")
  Matrix::writeMM(mat, file.path(mtx_dir, "matrix.mtx"))
  writeLines(as.character(colnames(mat)), file.path(mtx_dir, "barcodes.tsv"))
  genes <- as.character(rownames(mat))
  # scanpy.read_10x_mtx expects 10x v2 genes.tsv (id, symbol).
  utils::write.table(
    data.frame(genes, genes, stringsAsFactors = FALSE),
    file.path(mtx_dir, "genes.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  mtx_dir
}

pick_primary <- function(sample, requested = "scdblfinder") {
  ran <- as.character(sample$detector[as.character(sample$status) == "ran"])
  if (requested %in% ran) {
    return(requested)
  }
  hit <- DETECTORS[DETECTORS %in% ran]
  if (length(hit)) {
    return(hit[[1]])
  }
  NULL
}

call_to_is_doublet <- function(calls) {
  calls <- as.character(calls)
  out <- rep(NA, length(calls))
  out[calls == "doublet"] <- TRUE
  out[calls == "singlet"] <- FALSE
  as.logical(out)
}

attach_doublet_columns <- function(x, cells, sample, primary = "scdblfinder") {
  bc <- if (inherits(x, "Seurat")) colnames(x) else colnames(x)
  cells$barcode <- as.character(cells$barcode)
  idx <- match(as.character(bc), cells$barcode)
  for (name in DETECTORS) {
    sc <- paste0(name, "_score")
    cc <- paste0(name, "_call")
    if (sc %in% names(cells)) {
      vals <- cells[[sc]][idx]
      if (inherits(x, "Seurat")) {
        x[[sc]] <- vals
      } else {
        SummarizedExperiment::colData(x)[[sc]] <- vals
      }
    }
    if (cc %in% names(cells)) {
      vals <- as.character(cells[[cc]][idx])
      vals[is.na(vals)] <- ""
      if (inherits(x, "Seurat")) {
        x[[cc]] <- vals
      } else {
        SummarizedExperiment::colData(x)[[cc]] <- vals
      }
    }
  }
    used <- pick_primary(sample, primary)
    if (!is.null(used)) {
      score <- cells[[paste0(used, "_score")]][idx]
      is_d <- call_to_is_doublet(cells[[paste0(used, "_call")]][idx])
      if (inherits(x, "Seurat")) {
        x[["doublet_score"]] <- score
        x[["predicted_doublet"]] <- is_d
        x[["is_doublet"]] <- is_d
        x@misc$doublet_primary <- used
      } else {
        SummarizedExperiment::colData(x)[["doublet_score"]] <- score
        SummarizedExperiment::colData(x)[["predicted_doublet"]] <- is_d
        SummarizedExperiment::colData(x)[["is_doublet"]] <- is_d
        md <- S4Vectors::metadata(x)
        md$doublet_primary <- used
        S4Vectors::metadata(x) <- md
      }
    }
    x
}

annotate_doublets <- function(
    x,
    output_dir = NULL,
    fast = FALSE,
    n_jobs = -1L,
    random_state = 42L,
    primary = "scdblfinder",
    python = Sys.getenv("DOUBLET_PYTHON", "python3")
) {
  mat <- counts_from_object(x)
  if (is.null(output_dir)) {
    output_dir <- tempfile("doublet_")
  }
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mtx_dir <- file.path(output_dir, "sample")
  write_counts_mtx(mat, mtx_dir)
  py <- file.path(annotate_script_dir(), "run_doublet_rate.py")
  if (!file.exists(py)) {
    stop("cannot find run_doublet_rate.py next to annotate_object.R")
  }
  cmd <- c(
    py, mtx_dir, "--output-dir", output_dir,
    "--n-jobs", as.character(n_jobs),
    "--random-state", as.character(as.integer(random_state))
  )
  if (isTRUE(fast)) {
    cmd <- c(cmd, "--fast")
  }
  st <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) {
    stop(paste(c("run_doublet_rate.py failed:", st), collapse = "\n"))
  }
  cells_path <- file.path(output_dir, "sample.doublet_cells.tsv")
  sample_path <- file.path(output_dir, "sample.doublet_sample.tsv")
  if (!file.exists(cells_path)) {
    stop("missing ", cells_path)
  }
  cells <- utils::read.delim(cells_path, check.names = FALSE, stringsAsFactors = FALSE)
  sample <- utils::read.delim(sample_path, check.names = FALSE, stringsAsFactors = FALSE)
  attach_doublet_columns(x, cells, sample, primary = primary)
}
