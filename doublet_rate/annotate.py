"""Write Detector Scores/Calls back onto AnnData (Scanpy) objects."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io_counts import as_counts_adata


def calls_to_is_doublet(calls) -> pd.Series:
    """True/False from Call; empty Call is NA (not a singlet)."""
    s = pd.Series(calls, dtype=object).fillna("").astype(str)
    out = pd.Series(pd.NA, index=s.index, dtype="boolean")
    out.loc[s == "doublet"] = True
    out.loc[s == "singlet"] = False
    return out


def pick_primary(sample: pd.DataFrame, requested: str = "scdblfinder") -> str | None:
    from .run_doublet_rate import ALWAYS, GATED

    ran = []
    if "status" in sample.columns and "detector" in sample.columns:
        ran = sample.loc[sample["status"].astype(str) == "ran", "detector"].astype(str).tolist()
    if requested in ran:
        return requested
    for name in ALWAYS + GATED:
        if name in ran:
            return name
    return None


def attach_obs(adata, cell: pd.DataFrame, sample: pd.DataFrame, *, primary: str = "scdblfinder") -> str | None:
    """Copy per-Detector columns plus Scanpy aliases ``doublet_score`` / ``predicted_doublet``.

    Aliases are one Detector (default scDblFinder), not a fused consensus.
    ``predicted_doublet`` is the Scanpy boolean; ``is_doublet`` is the same column.
    """
    from .run_doublet_rate import ALWAYS, GATED

    mapped = cell.copy()
    mapped["barcode"] = mapped["barcode"].astype(str)
    mapped = mapped.drop_duplicates("barcode").set_index("barcode")
    mapped = mapped.reindex(adata.obs_names.astype(str))
    for name in ALWAYS + GATED:
        score_col = f"{name}_score"
        call_col = f"{name}_call"
        if score_col in mapped.columns:
            adata.obs[score_col] = mapped[score_col].to_numpy()
        if call_col in mapped.columns:
            adata.obs[call_col] = mapped[call_col].fillna("").astype(str).to_numpy()
    used = pick_primary(sample, primary)
    if used is None:
        return None
    adata.obs["doublet_score"] = adata.obs[f"{used}_score"]
    is_d = calls_to_is_doublet(adata.obs[f"{used}_call"])
    is_d.index = adata.obs.index
    adata.obs["predicted_doublet"] = is_d
    adata.obs["is_doublet"] = adata.obs["predicted_doublet"]
    adata.uns["doublet_primary"] = used
    return used


def detect_doublets(
    adata,
    *,
    output_dir: str | Path | None = None,
    fast: bool = False,
    n_jobs: int = -1,
    random_state: int = 42,
    primary: str = "scdblfinder",
    inplace: bool = True,
    stem: str = "sample",
):
    """Scanpy-native entry: score one AnnData Sample and write results into ``obs``.

    Returns the same AnnData (or a copy if ``inplace=False``). Reads raw counts
    from ``layers['counts']`` if present, else ``.X``. Does not ask the caller
    to extract a matrix, does not use ``obsm`` embeddings, and does not remove
    cells.

    Writes::

        adata.obs["doublet_score"]
        adata.obs["predicted_doublet"]

    ``random_state`` (default 42) is used for PCA, neighbor graphs, sampling,
    and classifiers in every Detector.
    """
    from .run_doublet_rate import score_adata

    target = adata if inplace else adata.copy()
    work = as_counts_adata(target)
    if output_dir is None:
        import tempfile

        out_dir = Path(tempfile.mkdtemp(prefix="doublet_"))
    else:
        out_dir = Path(output_dir).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
    cell, sample = score_adata(
        work, out_dir, stem=stem, fast=fast, n_jobs=n_jobs, random_state=random_state
    )
    attach_obs(target, cell, sample, primary=primary)
    target.uns["doublet_random_state"] = int(random_state)
    return target


annotate_anndata = detect_doublets
