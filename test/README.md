# test/

Bundled Sample for Quick Start: 10x Genomics **pbmc3k** filtered gene-barcode MTX (2,700 cells × 32,738 genes, already cell-called raw UMI counts).

```
test/
  matrix.mtx
  barcodes.tsv
  genes.tsv
```

```bash
python3 -m doublet_rate test --output-dir test_out
```

This runs the full roster (n ≤ 20,000). Solo trains scVI and is slow.
