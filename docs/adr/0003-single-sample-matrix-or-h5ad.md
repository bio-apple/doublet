# One Sample per run, as 10x MTX/H5 or h5ad raw counts

The skill accepts one Sample: a 10x MTX directory, a 10x `*.h5`, or a single-sample h5ad of already cell-called raw RNA counts. It does not split merged objects, does not run on integrated embeddings, does not read csv/tsv expression tables, and does not consume hashing, genotype, ATAC, or protein assays. Doublets cannot form across independently captured Samples, so a multi-sample object would be the wrong analysis unit.
