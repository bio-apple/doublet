# Architecture Decision Records

This directory is the **decision log** for the RNA-only doublet-rate tool. Each file records one choice that is hard to reverse, surprising without context, and the result of a real trade-off.

These are not tutorials and not a changelog. Read them when the code looks oddly strict — for example it reports a rate but never drops cells, or it refuses `nExp` / `dbr` as a Call cutoff — and you want the reason.

The glossary for the words used here is [CONTEXT.md](../../CONTEXT.md). How to run the tool is [SKILL.md](../../SKILL.md). Parameters and citations are [reference.md](../../reference.md).

## How to read an ADR

Each file is one paragraph: context, decision, why. Status is accepted unless a later ADR says otherwise. Numbering is sequential (`0001`, `0002`, …) and is not reused.

## Index

| ADR | Decision |
|---|---|
| [0001](0001-score-call-rate-not-removal.md) | Report Score, Call, and Predicted Doublet Rate. Do not remove cells. |
| [0002](0002-no-expected-rate-for-calls.md) | Do not use Expected Doublet Rate (`nExp`, `dbr`, 10x loading formula, Solo calibration) to place a Call. |
| [0003](0003-single-sample-matrix-or-h5ad.md) | One already cell-called RNA Sample per run (10x MTX / h5 / h5ad, or the same counts in AnnData, Seurat, or SingleCellExperiment). |
| [0004](0004-parallel-methods-no-fusion.md) | Report each Detector side by side. Do not fuse Calls into a consensus. |
| [0005](0005-detector-roster-and-size-gate.md) | Roster of score-capable methods; skip DoubletDecon; size-gate (and `--fast`-skip) DoubletDetection, DoubletFinder, and Solo. |
| [0006](0006-native-then-mad-calls.md) | Native data-driven threshold when it exists (Scrublet, scDblFinder `dbr.sd=1`); otherwise Griffiths/MAD. Never top-x% of Scores. |
| [0007](0007-no-qc-clustering-integration.md) | Do not QC, cluster, annotate, or integrate. The Sample is already cell-called raw counts. |
| [0008](0008-mixed-r-python.md) | Run R and Python Detectors in their native languages and join on barcode. Python `doublet_rate` orchestrates. |
| [0009](0009-continue-on-detector-failure.md) | One Detector failing does not abort the Sample; it is recorded as skipped. |
| [0010](0010-scripts-enforce-call-rules.md) | Do not reimplement Detector calls in a one-off notebook. Superseded in part by 0013: contract is the packages; `scripts/` are shims. |
| [0011](0011-rate-denominator-is-n-called.md) | Predicted Doublet Rate = `n_doublet / n_called`. Empty Calls are not singlets. |
| [0012](0012-primary-detector-scanpy-columns.md) | `doublet_score` / `predicted_doublet` / `is_doublet` copy one Primary Detector (default scDblFinder), not a consensus. |
| [0013](0013-dual-r-python-packages.md) | R package `doubletRate` and Python package `rna-doublet-rate` (`doublet_rate`) in one repo; `scripts/` are shims. |
