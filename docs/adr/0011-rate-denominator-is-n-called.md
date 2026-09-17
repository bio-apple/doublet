# Predicted Doublet Rate uses n_called, not n_input

The denominator is cells with a non-empty Call from that Detector. Treating missing Calls as singlets dilutes the rate when a Detector drops barcodes; treating them as doublets invents doublets. Record `n_input`, `n_called`, and `n_doublet` so the fraction is auditable.