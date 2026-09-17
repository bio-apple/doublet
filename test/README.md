# test/

Bundled Sample for Quick Start: 10x Genomics **pbmc3k** filtered gene-barcode MTX (2,700 cells × 32,738 genes, already cell-called raw UMI counts).

```
test/
  matrix.mtx
  barcodes.tsv
  genes.tsv
```

Smoke test (skips DoubletDetection, DoubletFinder, and Solo):

```bash
python3 -m doublet_rate test --output-dir test_out --fast
```

Those three rows in `test.doublet_sample.tsv` are `skipped_reason=fast:skip_gated`. Omit `--fast` for the full roster.