\name{annotate_doublets}
\alias{annotate_doublets}
\title{Annotate a Seurat or SingleCellExperiment Sample}
\usage{
annotate_doublets(
  x,
  output_dir = NULL,
  fast = FALSE,
  n_jobs = -1L,
  random_state = 42L,
  primary = "scdblfinder",
  python = Sys.getenv("DOUBLET_PYTHON", "python3")
)
}
\arguments{
\item{x}{A Seurat or SingleCellExperiment object with raw counts.}
\item{output_dir}{Optional directory for TSV side outputs. Default is a tempfile.}
\item{fast}{Skip DoubletDetection, DoubletFinder, and Solo
(\code{skipped_reason=fast:skip_gated}).}
\item{n_jobs}{Worker count (\code{-1} = all CPUs).}
\item{random_state}{Seed for PCA, neighbors, sampling, and classifiers.}
\item{primary}{Detector copied to \code{doublet_score} / \code{predicted_doublet} /
\code{is_doublet}. Not a consensus.}
\item{python}{Python executable that can run \code{python -m doublet_rate}
(default \code{DOUBLET_PYTHON} or \code{python3}).}
}
\value{
\code{x} with per-Detector \code{\{name\}_score} / \code{\{name\}_call} plus
convenience columns. Cells are not removed. The copied Detector name is stored
as Seurat \code{misc$doublet_primary} or SCE \code{metadata$doublet_primary}.
}
\description{
Runs the Detector roster by calling \code{python -m doublet_rate}. Requires the
Python package \code{doublet_rate}. Missing detector packages skip that Detector;
missing Python is a hard failure. \code{doublet_score}, \code{predicted_doublet},
and \code{is_doublet} copy one Primary Detector (default scDblFinder), not a
consensus.
}
