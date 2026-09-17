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


DEFAULT_RANDOM_STATE = 42


def seed_everything(random_state: int = DEFAULT_RANDOM_STATE) -> int:
    """Seed Python, NumPy, and Torch so PCA / neighbors / sampling / classifiers agree."""
    import random

    seed = int(random_state)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch

        torch.manual_seed(seed)
        cuda = getattr(torch, "cuda", None)
        if cuda is not None and cuda.is_available():
            cuda.manual_seed_all(seed)
    except ImportError:
        pass
    return seed


def resolve_n_jobs(n_jobs: int = -1) -> int:
    """Positive worker count. ``-1`` / ``None`` → all CPUs (sklearn/joblib convention)."""
    import os

    cpus = os.cpu_count() or 1
    if n_jobs is None or n_jobs < 0:
        return max(1, cpus)
    return max(1, int(n_jobs))


def knn_indices(
    X,
    k: int,
    *,
    metric: str = "euclidean",
    n_jobs: int = 1,
    random_state: int = DEFAULT_RANDOM_STATE,
) -> np.ndarray:
    """k neighbors per row, excluding self. pynndescent when present, else sklearn."""
    X = np.asarray(X)
    n = int(X.shape[0])
    k = max(1, min(int(k), n - 1))
    n_jobs = resolve_n_jobs(n_jobs)
    try:
        from pynndescent import NNDescent

        index = NNDescent(
            X,
            n_neighbors=min(k + 1, n),
            metric=metric,
            n_jobs=n_jobs,
            random_state=int(random_state),
            low_memory=True,
        )
        knn = np.asarray(index.neighbor_graph[0], dtype=int)
        return _knn_drop_self(knn, k)
    except Exception:
        from sklearn.neighbors import NearestNeighbors

        algo = "brute" if metric == "cosine" else "auto"
        nbrs = NearestNeighbors(
            n_neighbors=k + 1,
            metric=metric,
            algorithm=algo,
            n_jobs=n_jobs,
        ).fit(X)
        return _knn_drop_self(np.asarray(nbrs.kneighbors(return_distance=False), dtype=int), k)


def _knn_drop_self(knn: np.ndarray, k: int) -> np.ndarray:
    n = knn.shape[0]
    if knn.shape[1] > k and np.array_equal(knn[:, 0], np.arange(n)):
        return knn[:, 1 : k + 1]
    out = np.empty((n, k), dtype=int)
    for i in range(n):
        row = knn[i][knn[i] != i]
        out[i] = row[:k]
    return out


def _torch_accelerator() -> str | None:
    """CUDA (``gpu``) or Apple MPS. ``None`` if neither is present — Solo does not train on CPU."""
    try:
        import torch
    except ImportError:
        return None
    if torch.cuda.is_available():
        return "gpu"
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        return "mps"
    return None


def run_scrublet(adata, outdir: Path, random_state: int = DEFAULT_RANDOM_STATE, n_jobs: int = -1) -> None:
    name = "scrublet"
    n_jobs = resolve_n_jobs(n_jobs)
    random_state = seed_everything(random_state)
    try:
        import scrublet as scr
    except ImportError:
        _skip(outdir, name, "missing_package:scrublet")
        return
    try:
        counts = _counts_csr(adata)
        # expected_doublet_rate is a constructor default; it is not used to place the Call.
        import scrublet.helper_functions as hf

        orig_graph = hf.get_knn_graph

        def _graph(X, k=5, dist_metric="euclidean", approx=False, return_edges=True, random_seed=0):
            knn = knn_indices(
                X,
                k,
                metric=dist_metric,
                n_jobs=n_jobs,
                random_state=random_seed if random_seed is not None else random_state,
            )
            if not return_edges:
                return knn
            links = set()
            for i in range(knn.shape[0]):
                for j in knn[i]:
                    links.add(tuple(sorted((i, int(j)))))
            return links, knn

        hf.get_knn_graph = _graph
        try:
            scrub = scr.Scrublet(counts, random_state=random_state)
            scores, pred = scrub.scrub_doublets(verbose=False)
        finally:
            hf.get_knn_graph = orig_graph
        native = None
        if pred is not None:
            native = np.where(pred, "doublet", "singlet")
        _write(outdir, name, adata.obs_names, scores, native)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")


def run_doubletdetection(adata, outdir: Path, random_state: int = DEFAULT_RANDOM_STATE, n_jobs: int = -1) -> None:
    name = "doubletdetection"
    n_jobs = resolve_n_jobs(n_jobs)
    random_state = seed_everything(random_state)
    try:
        import doubletdetection
    except ImportError:
        _skip(outdir, name, "missing_package:doubletdetection")
        return
    try:
        counts = _counts_csr(adata)
        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            clustering_algorithm="leiden",
            standard_scaling=True,
            random_state=random_state,
            n_jobs=n_jobs,
            verbose=False,
        )
        clf.fit(counts)
        scores = np.asarray(clf.doublet_score(), dtype=float)
        _write(outdir, name, adata.obs_names, scores)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")


def run_solo(adata, outdir: Path, random_state: int = DEFAULT_RANDOM_STATE) -> None:
    name = "solo"
    try:
        import scvi
    except ImportError:
        _skip(outdir, name, "missing_package:scvi-tools")
        return
    accelerator = _torch_accelerator()
    if accelerator is None:
        _skip(outdir, name, "no_gpu")
        return
    try:
        from scvi.external import SOLO
        from scvi.model import SCVI

        ad = adata.copy()
        if not sparse.issparse(ad.X):
            ad.X = sparse.csr_matrix(np.asarray(ad.X))
        seed_everything(random_state)
        scvi.settings.seed = int(random_state)
        SCVI.setup_anndata(ad)
        vae = SCVI(ad)
        vae.train(accelerator=accelerator)
        solo = SOLO.from_scvi_model(vae)
        solo.train(accelerator=accelerator)
        pred = solo.predict(soft=True, include_simulated_doublets=False)
        if "doublet" not in pred.columns:
            _skip(outdir, name, "error:solo predict() has no 'doublet' column")
            return
        scores = pred["doublet"].to_numpy()
        barcodes = pred.index.astype(str) if getattr(pred, "index", None) is not None else ad.obs_names
        _write(outdir, name, barcodes, scores)
    except Exception as exc:
        _skip(outdir, name, f"error:{exc}")
