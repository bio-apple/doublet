# RNA-only doublet rate — reference

Companion to `SKILL.md`. The agent runs `scripts/run_doublet_rate.py`; this file records parameters, failure modes, and sources.

## Forbidden inputs

These place an Expected Doublet Rate into the Call and are not used:

| Parameter | Detector | Why it is forbidden |
|---|---|---|
| `nExp` as cutoff | DoubletFinder | Sets how many cells are labeled doublet; does not change pANN (Zhang et al. 2024) |
| `dbr` / `dbr.per1k` as cutoff | scDblFinder | 10x-style expected proportion |
| `expected_doublet_rate` as cutoff | Scrublet | Not used to place the Call here; Predicted Doublet Rate is not Scrublet's "overall" (detectable-fraction-adjusted) estimate |
| `--expected_number_of_doublets` / logit calibration | Solo | Bernstein et al. 2020 calibrate probabilities with a background rate |
| top x% of Scores | any | Writes a chosen identification rate into Predicted Doublet Rate (Xi and Li 2021 protocol) |
| rhop matching expected count | DoubletDecon | Out of roster (no continuous Score) |

scDblFinder still estimates a default `dbr` internally if omitted. `dbr.sd=1` disables that expectation at thresholding (scDblFinder docs): the Call is driven by misclassification of artificial doublets.

DoubletFinder is called with `nExp=1` only because the function requires the argument. The classification column is discarded.

## Call rules

**Native, no Expected Doublet Rate**

- Scrublet: `call_doublets(threshold=None)` / `scrub_doublets` automatic minimum between modes of simulated doublet scores (Wolock et al. 2019). Fallback: MAD if the automatic threshold is missing.
- scDblFinder: `scDblFinder(clusters=NULL, dbr.sd=1)` (Germain et al. 2021; OSCA.advanced ch. 8).

**MAD / Griffiths (everyone else)**

On that Detector's Score, `log1p`, median + 3 MAD (constant 1.4826), high tail = doublet. Same idea as OSCA converting doublet densities with `doubletThresholding(..., method="griffiths")` / Pijuan-Sala et al. 2019 large outliers. Applied to cxds, bcds, hybrid, DoubletFinder pANN, DoubletDetection scores, and Solo soft doublet scores.

Do not use DoubletDetection `predict(p_thresh=..., voter_thresh=...)` — those are global defaults, not a per-Sample score-distribution rule.

## Detector notes

| Detector | Ecosystem | Score | Notes |
|---|---|---|---|
| scDblFinder | R | `scDblFinder.score` | No clusters (OSCA allows `clusters=NULL`) |
| Scrublet | Python | `doublet_scores_obs_` | Constructor default expected rate is ignored for Calls |
| cxds | R / scds | `cxds_score` | Co-expression; no expected rate for scoring (Bais and Kostka 2020) |
| bcds | R / scds | `bcds_score` | |
| hybrid | R / scds | `hybrid_score` | Does not support specifying expected rate (Zhang et al. 2024) |
| DoubletDetection | Python | `BoostClassifier.doublet_score()` | Skip if n>20,000 (Xi and Li 2021: poor scaling) |
| DoubletFinder | R | pANN only | Internal PCA is method machinery, not a skill product. Skip if n>20,000 |
| Solo | Python / scvi-tools | soft `doublet` probability, uncalibrated | Train per Sample on CUDA or MPS if present, else CPU. Skip if n>20,000 (Demuxafy: median ~13 h at ~20k) |
| DoubletDecon | — | — | Not run: binary output, no Score (Xi and Li 2021) |

Size gate: `n_input > 20000`.

## I/O

**Accept**

- 10x MTX directory
- 10x `filtered_feature_bc_matrix.h5` (`gex_only=True`)
- h5ad with raw counts in `layers['counts']` or `.X`

**Reject**

- csv/tsv
- negative or clearly non-integer expression
- `obs` sample/batch columns with more than one value
- empty matrices

**Outputs**

Cell table: all input barcodes. Missing Detector values stay empty.

Sample table: one row per Detector in the roster, including skips.

`predicted_doublet_rate = n_doublet / n_called` (ADR 0011). Report `n_input`, `n_scored`, `n_called` so the fraction is auditable.

A Detector crash writes `{detector}.skip.txt` and continues (ADR 0009).

## Runtime

Mixed R + Python (ADR 0008). The driver exports gene-by-cell Matrix Market for R, runs `detectors_r.R` and the Python Detectors, then joins on barcode.

Seeds: Python `random_state=0`, R `set.seed(0)`.

## Failure modes

| Symptom | Likely cause | Handling |
|---|---|---|
| `matrix does not look like raw UMI counts` | log-normalized `.X` without `layers['counts']` | Fix the h5ad; do not run |
| `obs['batch'] has N values` | merged object | Split outside this skill |
| Scrublet `native_call` missing | unimodal simulated scores | MAD fallback |
| DoubletFinder skip on Seurat 5 layers | API mismatch | Record error; other Detectors still run |
| Solo skip | no GPU/CPU time, no scvi-tools | Expected; rate still comes from other Detectors |
| All R Detectors skip | `Rscript` missing or R packages missing | Install via `scripts/install_r_packages.R` |
| Trajectory-like intermediates called doublet | OSCA 8.5 | Interpret; do not auto-remove |

## Sources

- OSCA.advanced 3.23 ch. 8, [Doublet detection](https://bioconductor.org/books/3.23/OSCA.advanced/doublet-detection.html)
- Xi and Li, *Cell Systems* 2021 (©2020). Benchmark of computational doublet-detection methods
- Xi and Li, *STAR Protocols* 2021. Protocol for eight methods
- Germain et al., *F1000Research* 2021. scDblFinder
- Bais and Kostka, *Bioinformatics* 2020. scds (cxds/bcds/hybrid)
- Wolock et al., *Cell Systems* 2019. Scrublet
- McGinnis et al., *Cell Systems* 2019. DoubletFinder
- Gayoso et al. DoubletDetection
- Bernstein et al., *Cell Systems* 2020. Solo
- Neavin et al., *Genome Biology* 2024. Demuxafy (runtime at ~20k; genotype demux is out of scope)
- Zhang et al., *Cell Genomics* 2024. Expected rate changes DoubletFinder calls, not scores; hybrid has no expected-rate argument
- She et al., *CSBJ* 2025. Expected rate as a removal-count knob — not used here
- Liu et al., *Briefings in Bioinformatics* 2025. OmniDoublet (multimodal; RNA-only subset of comparisons only)

Local PDFs of the 2020–2025 papers live in this repository. Do not use hashing/genotype methods from OSCA 8.4 or Demuxafy, and do not use OmniDoublet's multimodal fusion.
