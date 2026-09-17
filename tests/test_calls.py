"""Checks for MAD Calls, Predicted Doublet Rate denominator, and n_jobs."""

from __future__ import annotations

import os
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
import pandas as pd

from doublet_rate.calls import griffiths_mad_calls, rate_from_calls
from doublet_rate.detectors_python import DEFAULT_RANDOM_STATE, resolve_n_jobs, seed_everything


def test_mad_calls_high_outliers():
    rng = np.random.default_rng(0)
    scores = np.concatenate([rng.normal(1.0, 0.1, 200), np.array([20.0, 22.0, 25.0])])
    calls = griffiths_mad_calls(scores)
    assert (calls[-3:] == "doublet").all()
    assert (calls[:-3] == "singlet").sum() >= 180


def test_empty_not_singlet():
    calls = pd.Series(["doublet", "singlet", "", "doublet"])
    out = rate_from_calls(calls, n_scored=3)
    assert out["n_called"] == 3
    assert out["n_doublet"] == 2
    assert abs(out["predicted_doublet_rate"] - 2 / 3) < 1e-12


def test_resolve_n_jobs():
    assert resolve_n_jobs(1) == 1
    assert resolve_n_jobs(4) == 4
    assert resolve_n_jobs(0) == 1
    n = resolve_n_jobs(-1)
    assert n == max(1, os.cpu_count() or 1)
    assert resolve_n_jobs(None) == n


def test_seed_everything_is_repeatable():
    assert DEFAULT_RANDOM_STATE == 42
    seed_everything(42)
    a = np.random.rand(4)
    seed_everything(42)
    b = np.random.rand(4)
    assert np.array_equal(a, b)
    seed_everything(1)
    c = np.random.rand(4)
    assert not np.array_equal(a, c)


if __name__ == "__main__":
    test_mad_calls_high_outliers()
    test_empty_not_singlet()
    test_resolve_n_jobs()
    test_seed_everything_is_repeatable()
    print("ok")
