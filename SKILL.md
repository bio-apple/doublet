---
name: rna-doublet-rate
description: Runs RNA-only multi-detector doublet scoring, data-driven calls, and predicted doublet rates on one 10x MTX/H5 or h5ad sample without removing cells. Use when the user asks for doublet detection, doublet rate, 多胞, scDblFinder, Scrublet, cxds, bcds, hybrid, Solo, DoubletDetection, or DoubletFinder on scRNA-seq.
---

# RNA-only doublet rate — run contract

This file is the **run contract**. CLI flags, output tables, the Detector roster, and object APIs are defined here, not in the README. The README is only the five-minute start plus a documentation index.

| Need | File |
|---|---|
| Five-minute start | [README.md](README.md) |
| Words this repo uses (Sample, Detector, Score, Call, rate) | [CONTEXT.md](CONTEXT.md) |
| Exact thresholds, skip reasons, citations | [reference.md](reference.md) |
| Why these rules exist | [docs/adr/](docs/adr/) |

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
python -m doublet_rate INPUT --n-jobs 8
python -m doublet_rate INPUT --random-state 42
python -m doublet_rate INPUT --primary scdblfinder
python -m doublet_rate INPUT --write-h5ad
# shim: python scripts/run_doublet_rate.py INPUT
# same CLI: doublet-rate INPUT
```

Install check (roster only, no Sample):

```bash
python -m doublet_rate --list-detectors
```

Full roster on [`test/`](test/) (Solo is slow):

```bash
python -m doublet_rate test --output-dir test_out --write-h5ad
```

| Flag | Meaning |
|---|---|
| `--list-detectors` | Print the roster (`always` / `gated`) and exit. No Sample. Does not probe whether packages are installed |
| `--output-dir` | Directory for TSV (and optional h5ad) outputs. Default: beside the input path; for a directory input, the parent of that directory |
| `--n-jobs` | Workers for artificial-doublet / KNN steps. Default `-1` = all CPUs. Passed to scDblFinder `BPPARAM`, DoubletFinder `paramSweep(num.cores)`, Scrublet KNN, DoubletDetection |
| `--random-state` | Seed for PCA, neighbor graphs, sampling, and classifiers. Default `42`. Same value in every Detector. Python: `detect_doublets(adata, random_state=42)`. R: `annotate_doublets(x, random_state = 42)` |
| `--write-h5ad` | Write `{stem}.doublet.h5ad` with `obs['doublet_score']`, `obs['predicted_doublet']`, `obs['is_doublet']`, plus per-Detector columns |
| `--primary` | Detector copied into those convenience columns (default `scdblfinder`). Not a consensus. If it did not run, copy the first Detector that did |

## Input

`INPUT` is one Sample of already cell-called raw RNA counts:

| | Required |
|---|---|
| Assay | RNA only (10x MTX directory, 10x `.h5`, or single-sample `.h5ad`) |
| Values | **Raw counts** (UMI/read integers). Use `layers['counts']` in h5ad if `.X` is normalized. Non-negative and ≥80% approximately integer, or the matrix is rejected |
| Barcodes | Unique cell barcodes |
| Unit | **One capture / one Sample**. Multi-sample `obs` columns (`sample`, `batch`, `orig.ident`, …) are rejected |
| Cell calling | **Already done**. Pass filtered barcodes, not empty droplets |

Protein, ATAC, hashing, and genotype assays are ignored. Merged objects, integrated embeddings, csv/tsv matrices, `obsm` embeddings, and log-normalized values are rejected.

This tool does **not** QC, filter low-quality cells, cluster, or integrate. Extra QC (mitochondrial fraction, min genes, …) is optional and must happen **before** you run this.

## Output

| File | Contents |
|---|---|
| `{stem}.doublet_cells.tsv` | `barcode` plus `{detector}_score` and `{detector}_call` |
| `{stem}.doublet_sample.tsv` | per Detector `n_input`, `n_scored`, `n_called`, `n_doublet`, `predicted_doublet_rate`, `status`, `skipped_reason`, `call_rule` |
| `{stem}.doublet.h5ad` | only with `--write-h5ad`; AnnData with the same columns on `obs` |

`status` is `ran`, `skipped` (size gate, missing package), `skipped:no_gpu` (Solo without CUDA or Apple MPS), or `failed` (crash, `skipped_reason=error:…`).

One Detector failing does not abort the Sample.

## In-memory objects

Do not convert the user through MTX by hand. Do not score `obsm` embeddings. Do not fuse Detectors into the convenience columns. Object APIs write TSV side output to `output_dir` when given, otherwise a tempfile.

R first (Seurat / SingleCellExperiment; `counts` assay / RNA counts layer). Requires `python -m doublet_rate`:

```r
library(doubletRate)
seu <- annotate_doublets(seu)
sce <- annotate_doublets(sce)
```

Python (AnnData; reads `layers['counts']` else `.X`):

```python
from doublet_rate import detect_doublets

adata = detect_doublets(adata)
```

`doublet_score`, `predicted_doublet`, and `is_doublet` copy **one** Detector (`primary`, default scDblFinder) on AnnData `obs`, Seurat `meta.data`, and SCE `colData`. Empty Calls are NA in `predicted_doublet`, not False. `is_doublet` is the same boolean as `predicted_doublet`. Per-Detector `{name}_score` / `{name}_call` are also written. The copied Detector name is `uns['doublet_primary']` (AnnData), `misc$doublet_primary` (Seurat), or `metadata$doublet_primary` (SCE). Cells are not removed.

## Detectors

| Detector | 默认（n≤20,000） | Size gate（n>20,000） |
|:---|:---:|:---:|
| scDblFinder | ✅ | ✅ |
| Scrublet | ✅ | ✅ |
| cxds | ✅ | ✅ |
| bcds | ✅ | ✅ |
| hybrid | ✅ | ✅ |
| DoubletDetection | ✅ | ⛔ |
| DoubletFinder (pANN) | ✅ | ⛔ |
| Solo | ✅ | ⛔ |

⛔ is `skipped_reason=size_gate:>20000`. DoubletDetection, DoubletFinder, and Solo are Gated Detectors; they stay on the roster. DoubletDecon is never run (no continuous Score).

```bash
python -m doublet_rate --list-detectors
```

prints the same eight names with `always` or `gated`.

Solo trains on CUDA or Apple MPS only. No GPU → `status=skipped:no_gpu` (`skipped_reason=no_gpu`), not a CPU fallback and not `failed`. A train/predict crash is `status=failed`.

Crashes are `status=failed` with `skipped_reason=error:…`. Expected absences are `skipped` (size gate, missing package) or `skipped:no_gpu`. Cell-table columns stay empty.

Python package `doublet_rate` orchestrates the roster. R Detectors still run in R (`detectors_r.R`) and Python Detectors in Python; results join on `barcode`.

## Calls and rate

- Scrublet: simulated-score bimodality (`call_rule=native:scrublet-bimodality`). If that fails, MAD.
- scDblFinder: `dbr.sd=1`, no user `dbr` (`call_rule=native:scdblfinder-dbr.sd=1`)
- All others: Griffiths/MAD high outliers on that Detector's Score (`call_rule=mad-griffiths`)
- **Predicted Doublet Rate = n_doublet / n_called**. Empty Calls are not singlets.

DoubletFinder `nExp=1` is API-only. Discard the DF class column. Call from MAD on pANN. Expected loading density is never a Call cutoff.

## Dependencies

Python: `pip install -e ".[full]"` (package name `rna-doublet-rate`; or `scripts/requirements.txt`). R: `R CMD INSTALL .` plus `Rscript scripts/install_r_packages.R`. The R package still needs the Python package. Missing *detector* packages skip that Detector; they do not stop the rest. Missing `python -m doublet_rate` stops the R function.

## Do not

- Remove cells
- Feed Expected Doublet Rate into a Call
- Cluster or integrate so that scDblFinder/DoubletFinder can use labels
- Fuse Detectors into a consensus Call
- Run on spatial, plate-based, or multi-sample objects
- Rewrite Detector invocation in a one-off notebook
