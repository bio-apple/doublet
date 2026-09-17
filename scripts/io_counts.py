"""Load one Sample of already cell-called raw RNA counts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

SAMPLE_KEYS = ("sample", "Sample", "samp", "batch", "orig.ident", "donor")
SIZE_GATE = 20000


class InputError(ValueError):
    pass


def sample_stem(path: Path) -> str:
    p = path.resolve()
    if p.is_dir():
        if p.name in {
            "filtered_feature_bc_matrix",
            "raw_feature_bc_matrix",
            "filtered_gene_bc_matrices",
        }:
            return p.parent.name
        return p.name
    return p.stem


def _as_array(X):
    if sparse.issparse(X):
        return X.data
    return np.asarray(X).ravel()


def assert_raw_counts(X) -> None:
    data = _as_array(X)
    data = data[np.isfinite(data)]
    if data.size == 0:
        raise InputError("count matrix is empty")
    if np.any(data < 0):
        raise InputError("negative values: not raw counts")
    int_frac = float(np.mean(np.abs(data - np.round(data)) < 1e-6))
    if int_frac < 0.8:
        raise InputError("matrix does not look like raw UMI counts (non-integer values)")


def assert_single_sample(adata) -> None:
    for key in SAMPLE_KEYS:
        if key in adata.obs.columns:
            n = adata.obs[key].nunique(dropna=True)
            if n > 1:
                raise InputError(
                    f"obs['{key}'] has {n} values; this skill accepts one Sample per run"
                )


def _use_counts_layer(adata):
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]
    return adata


def load_sample(path: str | Path):
    import scanpy as sc

    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise InputError(f"input not found: {path}")

    if path.is_dir():
        mtx = list(path.glob("matrix.mtx*"))
        if not mtx:
            raise InputError(f"no matrix.mtx in directory: {path}")
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    elif path.suffix == ".h5ad":
        adata = sc.read_h5ad(path)
        adata = _use_counts_layer(adata)
    elif path.suffix == ".h5":
        adata = sc.read_10x_h5(path, gex_only=True)
    else:
        raise InputError("input must be a 10x MTX directory, a 10x .h5, or a .h5ad")

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise InputError("empty matrix")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    assert_single_sample(adata)
    assert_raw_counts(adata.X)
    return adata


def export_mtx(adata, outdir: Path) -> Path:
    """Write gene-by-cell Matrix Market files for R detectors."""
    outdir.mkdir(parents=True, exist_ok=True)
    X = adata.X.T.tocsc() if sparse.issparse(adata.X) else sparse.csc_matrix(np.asarray(adata.X).T)
    mmwrite(outdir / "matrix.mtx", X)
    (outdir / "barcodes.tsv").write_text("\n".join(adata.obs_names.astype(str)) + "\n")
    genes = adata.var_names.astype(str)
    (outdir / "genes.tsv").write_text("\n".join(genes) + "\n")
    return outdir
