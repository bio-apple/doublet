"""AnnData attach helpers: counts-layer isolation and Scanpy obs aliases."""

from __future__ import annotations

import numpy as np
import pandas as pd

from annotate import attach_obs, calls_to_is_doublet, pick_primary
from io_counts import InputError, as_counts_adata


def test_as_counts_adata_uses_layer_not_x():
    import anndata as ad

    x_norm = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=float)
    counts = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
    adata = ad.AnnData(X=x_norm.copy())
    adata.obs_names = ["c1", "c2"]
    adata.var_names = ["g1", "g2"]
    adata.layers["counts"] = counts.copy()
    work = as_counts_adata(adata)
    assert np.array_equal(np.asarray(work.X), counts)
    assert np.array_equal(np.asarray(adata.X), x_norm)


def test_as_counts_adata_rejects_normalized_x():
    import anndata as ad

    adata = ad.AnnData(X=np.array([[0.1, 0.2], [0.3, 0.4]], dtype=float))
    adata.obs_names = ["c1", "c2"]
    adata.var_names = ["g1", "g2"]
    try:
        as_counts_adata(adata)
    except InputError as exc:
        assert "raw UMI" in str(exc)
        return
    raise AssertionError("expected InputError")


def test_is_doublet_empty_is_na():
    s = calls_to_is_doublet(["doublet", "singlet", "", "doublet"])
    assert bool(s.iloc[0]) is True
    assert bool(s.iloc[1]) is False
    assert pd.isna(s.iloc[2])


def test_attach_obs_aliases_primary_not_fusion():
    import anndata as ad

    adata = ad.AnnData(X=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]))
    adata.obs_names = ["c1", "c2", "c3"]
    adata.var_names = ["g1", "g2"]
    cell = pd.DataFrame(
        {
            "barcode": ["c1", "c2", "c3"],
            "scdblfinder_score": [0.1, 0.9, 0.2],
            "scdblfinder_call": ["singlet", "doublet", "singlet"],
            "scrublet_score": [0.4, 0.5, 0.6],
            "scrublet_call": ["singlet", "singlet", "doublet"],
        }
    )
    sample = pd.DataFrame(
        {
            "detector": ["scdblfinder", "scrublet", "solo"],
            "status": ["ran", "ran", "skipped"],
        }
    )
    used = attach_obs(adata, cell, sample, primary="scdblfinder")
    assert used == "scdblfinder"
    assert list(adata.obs["doublet_score"]) == [0.1, 0.9, 0.2]
    assert bool(adata.obs["predicted_doublet"].iloc[1]) is True
    assert bool(adata.obs["predicted_doublet"].iloc[0]) is False
    assert adata.obs["is_doublet"].equals(adata.obs["predicted_doublet"])
    assert adata.uns["doublet_primary"] == "scdblfinder"
    used2 = pick_primary(sample, "solo")
    assert used2 == "scdblfinder"


def test_h5ad_roundtrip_aliases():
    import tempfile
    from pathlib import Path

    import anndata as ad

    adata = ad.AnnData(X=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]))
    adata.obs_names = ["c1", "c2", "c3"]
    adata.var_names = ["g1", "g2"]
    cell = pd.DataFrame(
        {
            "barcode": ["c1", "c2", "c3"],
            "scdblfinder_score": [0.1, 0.9, 0.2],
            "scdblfinder_call": ["singlet", "doublet", ""],
        }
    )
    sample = pd.DataFrame({"detector": ["scdblfinder"], "status": ["ran"]})
    attach_obs(adata, cell, sample)
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "x.h5ad"
        adata.write_h5ad(path)
        b = ad.read_h5ad(path)
    assert list(np.asarray(b.obs["doublet_score"])) == [0.1, 0.9, 0.2]
    assert bool(b.obs["predicted_doublet"].iloc[1]) is True
    assert pd.isna(b.obs["predicted_doublet"].iloc[2])


if __name__ == "__main__":
    test_as_counts_adata_uses_layer_not_x()
    test_as_counts_adata_rejects_normalized_x()
    test_is_doublet_empty_is_na()
    test_attach_obs_aliases_primary_not_fusion()
    test_h5ad_roundtrip_aliases()
    print("ok")
