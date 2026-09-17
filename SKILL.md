---
name: rna-doublet-rate
description: Runs RNA-only multi-detector doublet scoring, data-driven calls, and predicted doublet rates on one 10x MTX/H5 or h5ad sample without removing cells. Use when the user asks for doublet detection, doublet rate, 多胞, scDblFinder, Scrublet, cxds, bcds, hybrid, Solo, DoubletDetection, or DoubletFinder on scRNA-seq.
---

# RNA-only doublet rate

Measure doublets on **one Sample**. Do not QC, cluster, integrate, fuse Detectors, or remove cells.

Execute `scripts/run_doublet_rate.py`. Do not reimplement Detector calls. Do not pass an Expected Doublet Rate as a Call cutoff (`nExp`, `dbr`, 10x 0.8%/1000 cells, Solo `--expected_number_of_doublets`, top-x% of scores).

## Run

```bash
python scripts/run_doublet_rate.py INPUT
python scripts/run_doublet_rate.py INPUT --output-dir DIR
```

`INPUT` is a 10x MTX directory (`matrix.mtx` + barcodes + features/genes), a 10x `*.h5`, or a single-sample `.h5ad` of already cell-called raw RNA counts. h5ad uses `layers['counts']` if present, else `.X`. Protein/ATAC/hashing/genotype assays are ignored. Merged objects, integrated embeddings, csv/tsv matrices, and log-normalized values are rejected.

Default output is beside the input:

- `{stem}.doublet_cells.tsv` — `barcode` plus `{detector}_score` and `{detector}_call`
- `{stem}.doublet_sample.tsv` — per Detector `n_input`, `n_scored`, `n_called`, `n_doublet`, `predicted_doublet_rate`, `status`, `skipped_reason`, `call_rule`

## Detectors

| Detector | ≤20,000 cells | >20,000 cells |
|---|---|---|
| scDblFinder, Scrublet, cxds, bcds, hybrid | run | run |
| DoubletDetection, DoubletFinder (pANN only), Solo | run | skip (`size_gate:>20000`) |
| DoubletDecon | never | never |

One Detector failing does not abort the Sample. Failed Detectors are `skipped` with a reason; their cell-table columns are empty.

R and Python run in their native languages and join on `barcode`.

## Calls and rate

- Scrublet: simulated-score bimodality (`call_rule=native:scrublet-bimodality`). If that fails, MAD.
- scDblFinder: `dbr.sd=1`, no user `dbr` (`call_rule=native:scdblfinder-dbr.sd=1`)
- All others: Griffiths/MAD high outliers on that Detector's Score (`call_rule=mad-griffiths`)
- **Predicted Doublet Rate = n_doublet / n_called**. Empty Calls are not singlets.

DoubletFinder `nExp=1` is API-only. Discard the DF class column. Call from MAD on pANN.

## Dependencies

Python: `scripts/requirements.txt`. R: `Rscript scripts/install_r_packages.R`. Missing packages skip that Detector; they do not stop the rest.

## Do not

- Remove cells
- Feed Expected Doublet Rate into a Call
- Cluster or integrate so that scDblFinder/DoubletFinder can use labels
- Fuse Detectors into a consensus Call
- Run on spatial, plate-based, or multi-sample objects
- Rewrite Detector invocation in a one-off notebook

Details: [reference.md](reference.md). Domain terms: [CONTEXT.md](CONTEXT.md).
