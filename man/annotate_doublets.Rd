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
\item{output_dir}{Optional directory for TSV side outputs.}
\item{fast}{Skip DoubletDetection, DoubletFinder, and Solo.}
\item{n_jobs}{Worker count (\code{-1} = all CPUs).}
\item{random_state}{Seed for PCA, neighbors, sampling, and classifiers.}
\item{primary}{Detector copied to \code{doublet_score} / \code{predicted_doublet}.}
\item{python}{Python executable.}
}
\value{
\code{x} with annotation columns added. Cells are not removed.
}
\description{
Runs the Detector roster. \code{doublet_score} and \code{predicted_doublet}
copy one Primary Detector (default scDblFinder), not a consensus.
}
