# Mixed R and Python, joined on barcode

The Detector list is split across ecosystems: R for scDblFinder, scds, and DoubletFinder; Python for Scrublet, Solo, and DoubletDetection. Run each in its native language and join results on cell barcode. Restricting to one language would delete methods; wrapping everything through a missing adapter would disguise failures as method choice. The orchestrator is the Python package `doublet_rate`; R `annotate_doublets` shells out to that CLI, which then calls `detectors_r.R` for the R Detectors.
