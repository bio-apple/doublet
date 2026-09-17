"""Checks for MAD Calls and Predicted Doublet Rate denominator."""

from __future__ import annotations

import numpy as np
import pandas as pd

from calls import griffiths_mad_calls, rate_from_calls


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


if __name__ == "__main__":
    test_mad_calls_high_outliers()
    test_empty_not_singlet()
    print("ok")
