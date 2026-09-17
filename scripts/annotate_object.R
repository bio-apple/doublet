# Shim: source() this file from a git checkout, or library(doubletRate) after install.
local({
  ofile <- NULL
  if (sys.nframe() >= 1) ofile <- sys.frame(1)$ofile
  root <- if (!is.null(ofile)) dirname(dirname(normalizePath(ofile))) else getwd()
  r_src <- file.path(root, "R", "annotate_doublets.R")
  if (!file.exists(r_src)) {
    stop("cannot find R/annotate_doublets.R; remotes::install_github('bio-apple/doublet')")
  }
  sys.source(r_src, envir = .GlobalEnv)
})
