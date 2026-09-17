"""Python Detectors: Scrublet, DoubletDetection, Solo. Scores only; native Call for Scrublet."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def _counts_csr(adata):
    X = adata.X
    if sparse.issparse(X):
        return X.tocsr()
    return sparse.csr_matrix(np.asarray(X))


def _write(outdir: Path, name: str, barcodes, score, native_call=None) -> None:
    df = pd.DataFrame({"barcode": np.asarray(barcodes, dtype=str), "score": np.asarray(score, dtype=float)})
    if native_call is not None:
        df["native_call"] = native_call
    df.to_csv(outdir / f"{name}.tsv", sep="\t", index=False)


def _skip(outdir: Path, name: str, reason: str) -> None:
    (outdir / f"{name}.skip.txt").write_text(reason + "\n")


def run_scrublet(adata, outdir: Path, seed: int = 0) -> None:
    name = "scrublet"
    try:
        import scrublet as scr
    except ImportError:
        _skip(outdir, name, "missing_package:scrublet")
        return
    try:
        counts = _counts_csr(adata)
        # expected_doublet_rate is a constructor default; it is not used to place the Call.
        scrub = scr.Scrublet(counts, random_state=seed)
        scores, pred = scrub.scrub_doublets(verbose=False)
        native = None
        if pred is not None:
            native = np.where(pred, "doublet", "singlet")
        _write(outdir, name, adata.obs_names, scores, native)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")


def run_doubletdetection(adata, outdir: Path, seed: int = 0) -> None:
    name = "doubletdetection"
    try:
        import doubletdetection
    except ImportError:
        _skip(outdir, name, "missing_package:doubletdetection")
        return
    try:
        counts = _counts_csr(adata)
        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            clustering_algorithm="louvain",
            standard_scaling=True,
            random_state=seed,
            n_jobs=1,
            verbose=False,
        )
        clf.fit(counts)
        scores = np.asarray(clf.doublet_score(), dtype=float)
        _write(outdir, name, adata.obs_names, scores)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")


def run_solo(adata, outdir: Path, seed: int = 0) -> None:
    name = "solo"
    try:
        import scvi
    except ImportError:
        _skip(outdir, name, "missing_package:scvi-tools")
        return
    try:
        from scvi.external import SOLO
        from scvi.model import SCVI

        ad = adata.copy()
        if not sparse.issparse(ad.X):
            ad.X = sparse.csr_matrix(np.asarray(ad.X))
        scvi.settings.seed = seed
        SCVI.setup_anndata(ad)
        vae = SCVI(ad)
        vae.train(accelerator="cpu")
        solo = SOLO.from_scvi_model(vae)
        solo.train(accelerator="cpu")
        pred = solo.predict(soft=True, include_simulated_doublets=False)
        if "doublet" in pred.columns:
            scores = pred["doublet"].to_numpy()
        else:
            scores = pred.iloc[:, -1].to_numpy()
        barcodes = pred.index.astype(str) if getattr(pred, "index", None) is not None else ad.obs_names
        _write(outdir, name, barcodes, scores)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")
