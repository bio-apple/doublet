# Dual installable packages in one repo; scripts/ stay shims

Seurat users install an R package (`doubletRate`); Scanpy users install a Python package (`rna-doublet-rate` / import `doublet_rate`). Both live in this repository so the Detector roster stays one contract. `scripts/` remains the agent entry (ADR 0010) as thin shims so existing commands keep working. Python Detectors are still invoked from R; missing Python packages become Skipped Detectors (ADR 0009), not a second roster.
