DETECTORS <- c(
  "scdblfinder", "scrublet", "cxds", "bcds", "hybrid",
  "doubletdetection", "doubletfinder", "solo"
)

repo_root <- function() {
  if (requireNamespace("doubletRate", quietly = TRUE)) {
    inst <- system.file(package = "doubletRate")
    parent <- dirname(inst)
    if (dir.exists(file.path(parent, "doublet_rate"))) {
      return(parent)
    }
  }
  ofile <- NULL
  n <- sys.nframe()
  if (n >= 1) {
    for (i in n:1) {
      ofile <- sys.frame(i)$ofile
      if (!is.null(ofile)) break
    }
  }
  if (!is.null(ofile)) {
    d <- dirname(normalizePath(ofile))
    if (basename(d) == "R") {
      return(dirname(d))
    }
    if (basename(d) == "scripts") {
      return(dirname(d))
    }
  }
  getwd()
}

python_cmd <- function() {
  Sys.getenv("DOUBLET_PYTHON", "python3")
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
  bc <- colnames(x)
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

run_python_cli <- function(mtx_dir, output_dir, n_jobs, random_state, python) {
  args <- c(
    "-m", "doublet_rate", mtx_dir,
    "--output-dir", output_dir,
    "--n-jobs", as.character(n_jobs),
    "--random-state", as.character(as.integer(random_state))
  )
  py_root <- repo_root()
  old <- Sys.getenv("PYTHONPATH", unset = NA_character_)
  if (dir.exists(file.path(py_root, "doublet_rate"))) {
    sep <- .Platform$path.sep
    cur <- if (is.na(old)) "" else old
    Sys.setenv(PYTHONPATH = if (nzchar(cur)) paste(py_root, cur, sep = sep) else py_root)
    on.exit({
      if (is.na(old)) Sys.unsetenv("PYTHONPATH") else Sys.setenv(PYTHONPATH = old)
    }, add = TRUE)
  }
  system2(python, args, stdout = TRUE, stderr = TRUE)
}

#' Annotate a Seurat or SingleCellExperiment Sample
#'
#' Runs the Detector roster via \code{python -m doublet_rate} and writes
#' per-Detector scores/calls plus \code{doublet_score} / \code{predicted_doublet}
#' / \code{is_doublet} from the Primary Detector. Does not remove cells.
#' Missing detector packages skip that Detector; missing Python is a hard failure.
#'
#' @param x A Seurat or SingleCellExperiment object with raw counts.
#' @param output_dir Optional directory for TSV side outputs. Default is a tempfile.
#' @param n_jobs Worker count (\code{-1} = all CPUs).
#' @param random_state Seed for PCA, neighbors, sampling, and classifiers.
#' @param primary Detector copied to \code{doublet_score} / \code{predicted_doublet} /
#'   \code{is_doublet}. Not a consensus.
#' @param python Python executable that can run \code{python -m doublet_rate}
#'   (default \code{DOUBLET_PYTHON} or \code{python3}).
#' @return \code{x} with annotation columns added. Cells are not removed.
#' @export
annotate_doublets <- function(
    x,
    output_dir = NULL,
    n_jobs = -1L,
    random_state = 42L,
    primary = "scdblfinder",
    python = python_cmd()
) {
  mat <- counts_from_object(x)
  if (is.null(output_dir)) {
    output_dir <- tempfile("doublet_")
  }
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mtx_dir <- file.path(output_dir, "sample")
  write_counts_mtx(mat, mtx_dir)
  st <- run_python_cli(mtx_dir, output_dir, n_jobs, random_state, python)
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) {
    stop(paste(c("python -m doublet_rate failed:", st), collapse = "\n"))
  }
  cells_path <- file.path(output_dir, "sample.doublet_cells.tsv")
  sample_path <- file.path(output_dir, "sample.doublet_sample.tsv")
  if (!file.exists(cells_path)) {
    stop("missing ", cells_path, "\n", paste(st, collapse = "\n"))
  }
  cells <- utils::read.delim(cells_path, check.names = FALSE, stringsAsFactors = FALSE)
  sample <- utils::read.delim(sample_path, check.names = FALSE, stringsAsFactors = FALSE)
  attach_doublet_columns(x, cells, sample, primary = primary)
}
