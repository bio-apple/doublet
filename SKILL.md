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

The contract is the installable packages: `python -m doublet_rate` / `doublet-rate`, `from doublet_rate import detect_doublets`, and `library(doubletRate); annotate_doublets(x)`. `scripts/run_doublet_rate.py` and `scripts/annotate_object.R` are shims. Do not reimplement Detector calls. The R function always shells out to `python -m doublet_rate`; missing Python is a hard failure, not a partial R-only run. Do not pass an Expected Doublet Rate as a Call cutoff (`nExp`, `dbr`, 10x 0.8%/1000 cells, Solo `--expected_number_of_doublets`, top-x% of scores).

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
python -m doublet_rate --list-detectors
python -m doublet_rate INPUT
python -m doublet_rate INPUT --output-dir DIR
python -m doublet_rate INPUT --fast
python -m doublet_rate INPUT --n-jobs 8
python -m doublet_rate INPUT --random-state 42
python -m doublet_rate INPUT --primary scdblfinder
python -m doublet_rate INPUT --write-h5ad
# shim: python scripts/run_doublet_rate.py INPUT
```

Install smoke test (skips the three Gated Detectors):

```bash
python -m doublet_rate test --output-dir test_out --fast
```

| Flag | Meaning |
|---|---|
| `--list-detectors` | Print the roster (`always` / `gated`) and exit. No Sample. Not an install check |
| `--fast` | Skip DoubletDetection, DoubletFinder, and Solo. Sample table: `status=skipped`, `skipped_reason=fast:skip_gated`. Install smoke test, not the size gate |
| `--n-jobs` | Workers for artificial-doublet / KNN steps. Default `-1` = all CPUs. Passed to scDblFinder `BPPARAM`, DoubletFinder `paramSweep(num.cores)`, Scrublet KNN, DoubletDetection |
| `--random-state` | Seed for PCA, neighbor graphs, sampling, and classifiers. Default `42`. Same value in every Detector |
| `--write-h5ad` | Write `{stem}.doublet.h5ad` with `obs['doublet_score']`, `obs['predicted_doublet']`, `obs['is_doublet']`, plus per-Detector columns |
| `--primary` | Detector copied into those convenience columns (default `scdblfinder`). Not a consensus. If it did not run, copy the first Detector that did |

## Input

`INPUT` is one Sample of already cell-called raw RNA counts:

- 10x MTX directory (`matrix.mtx` + barcodes + features/genes)
- 10x `*.h5`
- single-sample `.h5ad` (`layers['counts']` if present, else `.X`)

Counts must be non-negative and look like integers (≥80% of finite values within `1e-6` of an integer). Barcodes must be unique. Protein, ATAC, hashing, and genotype assays are ignored. Merged objects, integrated embeddings, csv/tsv matrices, `obsm` embeddings, and log-normalized values are rejected.

## Output

Default location is `--output-dir`. If omitted, files are written beside the input path (for a directory input: the parent of that directory).

| File | Contents |
|---|---|
| `{stem}.doublet_cells.tsv` | `barcode` plus `{detector}_score` and `{detector}_call` |
| `{stem}.doublet_sample.tsv` | per Detector `n_input`, `n_scored`, `n_called`, `n_doublet`, `predicted_doublet_rate`, `status`, `skipped_reason`, `call_rule` |
| `{stem}.doublet.h5ad` | only with `--write-h5ad`; AnnData with the same columns on `obs` |

## In-memory objects

Do not convert the user through MTX by hand. Do not score `obsm` embeddings. Do not fuse Detectors into the convenience columns. Object APIs write TSV side output to `output_dir` when given, otherwise a tempfile.

R first (Seurat / SingleCellExperiment; `counts` assay / RNA counts layer). Requires `python -m doublet_rate`:

```r
library(doubletRate)
x <- annotate_doublets(x, fast = TRUE, n_jobs = 8)
```

Python (AnnData; reads `layers['counts']` else `.X`):

```python
from doublet_rate import detect_doublets

adata = detect_doublets(adata)
```

`doublet_score`, `predicted_doublet`, and `is_doublet` copy **one** Detector (`primary`, default scDblFinder) on AnnData `obs`, Seurat `meta.data`, and SCE `colData`. Empty Calls are NA in `predicted_doublet`, not False. `is_doublet` is the same boolean as `predicted_doublet`. Cells are not removed.

## Detectors

| Detector | 默认（n≤20,000） | Size gate（n>20,000） | `--fast` |
|:---|:---:|:---:|:---:|
| scDblFinder | ✅ | ✅ | ✅ |
| Scrublet | ✅ | ✅ | ✅ |
| cxds | ✅ | ✅ | ✅ |
| bcds | ✅ | ✅ | ✅ |
| hybrid | ✅ | ✅ | ✅ |
| DoubletDetection | ✅ | ⛔ | ⛔ |
| DoubletFinder (pANN) | ✅ | ⛔ | ⛔ |
| Solo | ✅ | ⛔ | ⛔ |

⛔ in the Size gate column is `skipped_reason=size_gate:>20000`. ⛔ in `--fast` is `fast:skip_gated`. DoubletDetection, DoubletFinder, and Solo are Gated Detectors; they stay on the roster. DoubletDecon is never run (no continuous Score).

One Detector failing does not abort the Sample. Failed Detectors are `skipped` with a reason; their cell-table columns are empty.

Python package `doublet_rate` orchestrates the roster. R Detectors still run in R (`detectors_r.R`) and Python Detectors in Python; results join on `barcode`.

## Calls and rate

- Scrublet: simulated-score bimodality (`call_rule=native:scrublet-bimodality`). If that fails, MAD.
- scDblFinder: `dbr.sd=1`, no user `dbr` (`call_rule=native:scdblfinder-dbr.sd=1`)
- All others: Griffiths/MAD high outliers on that Detector's Score (`call_rule=mad-griffiths`)
- **Predicted Doublet Rate = n_doublet / n_called**. Empty Calls are not singlets.

DoubletFinder `nExp=1` is API-only. Discard the DF class column. Call from MAD on pANN.

## Dependencies

Python: `pip install -e ".[full]"` (package name `rna-doublet-rate`; or `scripts/requirements.txt`). R: `R CMD INSTALL .` plus `Rscript scripts/install_r_packages.R`. The R package still needs the Python package. Missing *detector* packages skip that Detector; they do not stop the rest. Missing `python -m doublet_rate` stops the R function.

## Do not

- Remove cells
- Feed Expected Doublet Rate into a Call
- Cluster or integrate so that scDblFinder/DoubletFinder can use labels
- Fuse Detectors into a consensus Call
- Treat `--fast` skips as a size-gate failure
- Run on spatial, plate-based, or multi-sample objects
- Rewrite Detector invocation in a one-off notebook
