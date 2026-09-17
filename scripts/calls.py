"""Data-driven Calls and Predicted Doublet Rate (no Expected Doublet Rate)."""

from __future__ import annotations

import numpy as np
import pandas as pd

NMADS = 3.0
MAD_CONSTANT = 1.4826


def griffiths_mad_threshold(scores: np.ndarray, nmads: float = NMADS) -> float:
    """High-outlier threshold on log1p(scores), matching scater::isOutlier(type='higher')."""
    x = np.asarray(scores, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("inf")
    x = np.log1p(np.clip(x, 0, None))
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)) * MAD_CONSTANT)
    if mad == 0:
        return float("inf")
    return med + nmads * mad


def griffiths_mad_calls(scores: np.ndarray, nmads: float = NMADS) -> np.ndarray:
    """Return singlet/doublet labels; non-finite scores stay empty."""
    raw = np.asarray(scores, dtype=float)
    out = np.full(raw.shape, "", dtype=object)
    finite = np.isfinite(raw)
    if not finite.any():
        return out
    x = np.log1p(np.clip(raw[finite], 0, None))
    thresh = griffiths_mad_threshold(raw, nmads=nmads)
    out[finite] = np.where(x > thresh, "doublet", "singlet")
    return out


def rate_from_calls(calls: pd.Series, n_scored: int) -> dict:
    """Predicted Doublet Rate = n_doublet / n_called. Empty Calls are not singlets."""
    called = calls.isin(["singlet", "doublet"])
    n_called = int(called.sum())
    n_doublet = int((calls == "doublet").sum())
    rate = (n_doublet / n_called) if n_called else float("nan")
    return {
        "n_scored": int(n_scored),
        "n_called": n_called,
        "n_doublet": n_doublet,
        "predicted_doublet_rate": rate,
    }
