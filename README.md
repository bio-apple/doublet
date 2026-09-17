# RNA-only doublet rate

Score every cell in **one Sample**, make a data-driven Call, and report a Predicted Doublet Rate (`n_doublet / n_called`). This repository **does not remove cells**.

## Input (read this first)

| | Required |
|---|---|
| Assay | RNA only (10x MTX directory, 10x `.h5`, or single-sample `.h5ad`) |
| Values | **Raw counts** (UMI/read integers). Use `layers['counts']` in h5ad if `.X` is normalized |
| Unit | **One capture / one Sample**. Multi-sample `obs` columns (`sample`, `batch`, `orig.ident`, …) are rejected |
| Cell calling | **Already done**. Pass filtered barcodes, not empty droplets |

This tool does **not** QC, filter low-quality cells, cluster, or integrate. Do cell calling upstream. Extra QC (mitochondrial fraction, min genes, …) is optional and must happen **before** you run this; we will not do it for you.

**Do not pass** log-normalized, scaled, or integrated matrices. They are rejected.

## Quick Start

```bash
git clone https://github.com/bio-apple/doublet.git
cd doublet
python3 -m pip install -r scripts/requirements.txt
Rscript scripts/install_r_packages.R
python3 examples/demo.py
```

`examples/demo.py` loads [`test/`](test/) (10x pbmc3k filtered MTX) and runs the fast Detectors (`--fast`). It should print `demo ok`.

Equivalent one-liner:

```bash
python3 scripts/run_doublet_rate.py test --output-dir test_out --fast
```

## Test data

[`test/`](test/) is the public 10x Genomics pbmc3k **filtered** gene-barcode matrix (2,700 cells × 32,738 genes, integer UMIs, already cell-called). Use it to verify the install and as the default example Sample.

Full roster including Solo (omit `--fast`) is slow: scVI/SOLO train for hundreds of epochs.

## Full Sample run

```bash
python3 scripts/run_doublet_rate.py PATH/TO/SAMPLE --output-dir OUT
```

Omit `--fast` to also run DoubletDetection, DoubletFinder (pANN), and Solo when the Sample has ≤ 20,000 cells. Solo uses CUDA or Apple MPS when present.

Outputs next to `--output-dir` (default: beside the input):

- `{stem}.doublet_cells.tsv` — per barcode, per Detector Score and Call
- `{stem}.doublet_sample.tsv` — per Detector `n_input`, `n_scored`, `n_called`, `n_doublet`, `predicted_doublet_rate`, `status`

One Detector failing does not abort the Sample; it is recorded as `skipped`.

## Detectors

| Detector | ≤20,000 cells | >20,000 cells |
|---|---|---|
| scDblFinder, Scrublet, cxds, bcds, hybrid | run | run |
| DoubletDetection, DoubletFinder (pANN), Solo | run | skip |
| `--fast` | skip the three slow Detectors | skip them |

Calls: Scrublet and scDblFinder use their native data-driven thresholds (no Expected Doublet Rate). Everyone else uses Griffiths/MAD high outliers. Expected loading density is never a Call cutoff.

## Docs

- Domain terms: [CONTEXT.md](CONTEXT.md)
- Agent/run contract: [SKILL.md](SKILL.md)
- Parameters and sources: [reference.md](reference.md)
