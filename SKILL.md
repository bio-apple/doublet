---
name: rna-doublet-rate
description: Runs RNA-only multi-detector doublet scoring, data-driven calls, and predicted doublet rates on one 10x MTX/H5 or h5ad sample without removing cells. Use when the user asks for doublet detection, doublet rate, 多胞, scDblFinder, Scrublet, cxds, bcds, hybrid, Solo, DoubletDetection, or DoubletFinder on scRNA-seq.
---

# RNA-only doublet rate — run contract

This file is the **run contract**. It tells a person or an agent what to execute, which inputs are legal, and what this tool must not do. It is not a glossary and not a literature review.

| Need | File |
|---|---|
| Words this repo uses (Sample, Detector, Score, Call, rate) | [CONTEXT.md](CONTEXT.md) |
| Parameters, I/O, failure modes, citations | [reference.md](reference.md) |
| Why these rules exist | [docs/adr/](docs/adr/) |
| Quick Start | [README.md](README.md) |

Measure doublets on **one Sample**. Do not QC, cluster, integrate, fuse Detectors, or remove cells.

Execute `scripts/run_doublet_rate.py` (or the AnnData / Seurat / SCE wrappers below). Do not reimplement Detector calls. Do not pass an Expected Doublet Rate as a Call cutoff (`nExp`, `dbr`, 10x 0.8%/1000 cells, Solo `--expected_number_of_doublets`, top-x% of scores).

## Contents

1. [Run](#run)
2. [Input](#input)
3. [Output](#output)
4. [In-memory objects](#in-memory-objects)
5. [Detectors](#detectors)
6. [Calls and rate](#calls-and-rate)
7. [Dependencies](#dependencies)
8. [Do not](#do-not)

## Run

```bash
python scripts/run_doublet_rate.py INPUT
python scripts/run_doublet_rate.py INPUT --output-dir DIR
python scripts/run_doublet_rate.py INPUT --fast
python scripts/run_doublet_rate.py INPUT --n-jobs 8
python scripts/run_doublet_rate.py INPUT --random-state 42
python scripts/run_doublet_rate.py INPUT --write-h5ad
```

| Flag | Meaning |
|---|---|
| `--fast` | Skip DoubletDetection, DoubletFinder, and Solo (install smoke test) |
| `--n-jobs` | Workers for artificial-doublet / KNN steps. Default `-1` = all CPUs. Passed to scDblFinder `BPPARAM`, DoubletFinder `paramSweep(num.cores)`, Scrublet KNN, DoubletDetection |
| `--random-state` | Seed for PCA, neighbor graphs, sampling, and classifiers. Default `42`. Same value in every Detector. |
| `--write-h5ad` | Write `{stem}.doublet.h5ad` with `obs['doublet_score']` and `obs['predicted_doublet']` plus per-Detector columns |
| `--primary` | Detector copied into those two Scanpy columns (default `scdblfinder`). Not a consensus |

`examples/demo.py` is the one-click check on [`test/`](test/).

## Input

`INPUT` is one Sample of already cell-called raw RNA counts:

- 10x MTX directory (`matrix.mtx` + barcodes + features/genes)
- 10x `*.h5`
- single-sample `.h5ad` (`layers['counts']` if present, else `.X`)

Protein, ATAC, hashing, and genotype assays are ignored. Merged objects, integrated embeddings, csv/tsv matrices, `obsm` embeddings, and log-normalized values are rejected.

## Output

Default location is beside the input (or `--output-dir`):

| File | Contents |
|---|---|
| `{stem}.doublet_cells.tsv` | `barcode` plus `{detector}_score` and `{detector}_call` |
| `{stem}.doublet_sample.tsv` | per Detector `n_input`, `n_scored`, `n_called`, `n_doublet`, `predicted_doublet_rate`, `status`, `skipped_reason`, `call_rule` |
| `{stem}.doublet.h5ad` | only with `--write-h5ad`; AnnData with the same columns on `obs` |

## In-memory objects

Do not convert the user through MTX by hand. Do not score `obsm` embeddings. Do not fuse Detectors into the two Scanpy columns.

Python (AnnData; reads `layers['counts']` else `.X`; does not replace `.X` when a counts layer exists):

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts").resolve()))
from annotate import detect_doublets

adata = detect_doublets(adata)
```

R (Seurat or SingleCellExperiment; `counts` assay / RNA counts layer):

```r
source("scripts/annotate_object.R")
x <- annotate_doublets(x, fast = TRUE, n_jobs = 8)
```

`obs['doublet_score']` and `obs['predicted_doublet']` (and the same names on Seurat `meta.data` / SCE `colData`) copy **one** Detector (`primary`, default scDblFinder). Empty Calls are NA in `predicted_doublet`, not False. Cells are not removed. `is_doublet` is the same boolean as `predicted_doublet`.

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
