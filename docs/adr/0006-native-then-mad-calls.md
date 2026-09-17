# Native auto-threshold first, then MAD/Griffiths

When a Detector can threshold from the data without an Expected Doublet Rate, use that rule (Scrublet simulated-score bimodality; scDblFinder with `dbr.sd=1`). Otherwise call high outliers on that Detector's Score with the Griffiths/MAD rule OSCA uses for doublet densities. Do not call the top x% of Scores: that writes the Expected Doublet Rate back into the Predicted Doublet Rate.
