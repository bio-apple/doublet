# RNA-only Doublet Rate Analysis

This context is computational, RNA-only analysis of doublets in scRNA-seq: score every cell, make a data-driven call, and report a predicted rate. It does not remove cells.

## Language

**Doublet**:
An artifactual library generated from two cells captured together.
_Avoid_: Multiplet (that includes three or more cells)

**Sample**:
One droplet capture of already cell-called RNA counts, provided as a 10x MTX directory, a 10x `*.h5`, or a single-sample h5ad. The analysis unit is exactly one Sample.
_Avoid_: Merged object, integrated object, batch as a synonym, empty droplets, csv/tsv expression tables

**Detector**:
One computational method that emits a Doublet Score for every cell in the Sample. Detectors are reported side by side and are not fused.
_Avoid_: Ensemble, consensus method, pipeline

**Doublet Score**:
A per-cell continuous value from one Detector, higher meaning more doublet-like.
_Avoid_: Probability (unless the method actually emits a calibrated probability), rank

**Doublet Call**:
A per-cell singlet-or-doublet label obtained by a data-driven threshold on that Detector's Doublet Score. Expected Doublet Rate is not an input to this threshold.
_Avoid_: Filter, removal, classification as a synonym for the whole analysis

**Predicted Doublet Rate**:
The proportion of cells with a non-empty Call from one Detector that were called doublet (`n_doublet / n_called`). Empty Calls are not singlets. Each Detector reports its own rate; rates are not fused.
_Avoid_: Doublet rate (unqualified), identification rate, expected rate, consensus rate

**Expected Doublet Rate**:
An a priori proportion taken from loading density or chemistry, supplied as an algorithm input. This context does not use it to make a Call.
_Avoid_: Using this term for Predicted Doublet Rate

**Skipped Detector**:
A Detector absent from this Sample's results because of the cell-count gate, because it cannot emit a Score under the rules of this context, or because it errored. Other Detectors still report.
_Avoid_: Failed QC, filtered method, aborted analysis

**Removal**:
Dropping called doublets from the dataset before downstream analysis. This context does not do Removal.
_Avoid_: Filtering doublets as if it were the analysis itself
