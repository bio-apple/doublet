# Bioconductor / CRAN packages for the R Detectors.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("SingleCellExperiment", "scDblFinder", "scds", "BiocNeighbors", "BiocParallel"), ask = FALSE, update = FALSE)
install.packages(c("Seurat", "Matrix"), repos = "https://cloud.r-project.org")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never")
