# RNA-only Doublet Rate Analysis

This file is the **glossary** for this repository. It defines the words used in the README, the run contract, and the decision log. It is not a user guide, not a command list, and not a place for software parameters.

| Need | File |
|---|---|
| Five-minute start | [README.md](README.md) |
| How to run (CLI, output, Detector roster) | [SKILL.md](SKILL.md) |
| Parameters, I/O, failures, citations | [reference.md](reference.md) |
| Why a rule exists | [docs/adr/](docs/adr/) |

## How the terms connect

A [Sample](#sample) is one capture of already cell-called RNA counts. Each [Detector](#detector) emits a [Doublet Score](#doublet-score) for every cell. A data-driven [Call Rule](#call-rule) turns that Score into a [Doublet Call](#doublet-call). The [Predicted Doublet Rate](#predicted-doublet-rate) is `n_doublet / n_called` for that Detector. Detectors are not fused. The [Primary Detector](#primary-detector) is the one Detector whose Score and Call are copied to the convenience columns; that copy is not a consensus. [Gated Detectors](#gated-detector) may be a [Skipped Detector](#skipped-detector) on a given Sample. [Expected Doublet Rate](#expected-doublet-rate) is not used to make a Call. [Removal](#removal) of cells is out of scope.

## Language

### Doublet

An artifactual library generated from two cells captured together.
_Avoid_: Multiplet (that includes three or more cells)

### Sample

One droplet capture of already cell-called RNA counts. On disk: a 10x MTX directory, a 10x `*.h5`, or a single-sample h5ad. In memory: AnnData, Seurat, or SingleCellExperiment carrying that same count matrix.
_Avoid_: Merged object, integrated object, batch as a synonym, empty droplets, csv/tsv expression tables

### Detector

One computational method that emits a Doublet Score for every cell in the Sample. Detectors are reported side by side and are not fused.
_Avoid_: Ensemble, consensus method, pipeline

### Doublet Score

A per-cell continuous value from one Detector, higher meaning more doublet-like.
_Avoid_: Probability (unless the method actually emits a calibrated probability), rank

### Doublet Call

A per-cell label obtained by a data-driven Call Rule on that Detector's Doublet Score: doublet, singlet, or empty. An empty Call is a third state, not a singlet. Expected Doublet Rate is not an input to this threshold.
_Avoid_: Filter, removal, classification as a synonym for the whole analysis

### Call Rule

The named data-driven procedure that turns one Detector's Score into a Call. Native when the method has a data-driven threshold that does not take Expected Doublet Rate; otherwise Griffiths/MAD. Each Detector has its own Call Rule. Call Rules are not fused.
_Avoid_: Cutoff, expected-rate threshold, consensus rule

### Predicted Doublet Rate

The proportion of cells with a non-empty Call from one Detector that were called doublet (`n_doublet / n_called`). Empty Calls are not singlets. Each Detector reports its own rate; rates are not fused.
_Avoid_: Doublet rate (unqualified), identification rate, expected rate, consensus rate

### Expected Doublet Rate

An a priori proportion taken from loading density or chemistry, supplied as an algorithm input. This context does not use it to make a Call.
_Avoid_: Using this term for Predicted Doublet Rate

### Gated Detector

DoubletDetection, DoubletFinder, and Solo. Members of the roster that become a Skipped Detector when the Sample is above the size gate.
_Avoid_: Optional detector as if they were outside the roster, failed method

### Skipped Detector

A Detector absent from this Sample's results because of the size gate, because Solo has no CUDA or Apple MPS, or because it cannot emit a Score under the rules of this context (including a missing package). Other Detectors still report.
_Avoid_: Failed QC, filtered method, aborted analysis, using this term for a crash

### Failed Detector

A Detector that was attempted and crashed. Other Detectors still report. A missing GPU for Solo is a Skipped Detector, not a Failed Detector.
_Avoid_: Skipped Detector as a synonym for a crash

### Primary Detector

The one Detector whose Score and Call are copied onto the convenience columns `doublet_score`, `predicted_doublet`, and `is_doublet` on AnnData, Seurat, and SingleCellExperiment. It is a member of the roster, not a fused Call. If the requested Primary Detector did not run, the first Detector that did run is copied instead.
_Avoid_: Consensus call, default method as if it were the only Detector, main label, Scanpy-only columns

### Removal

Dropping called doublets from the dataset before downstream analysis. This context does not do Removal.
_Avoid_: Filtering doublets as if it were the analysis itself
