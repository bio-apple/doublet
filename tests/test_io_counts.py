"""Checks for 1-based / OOB sparse index sanitization used by MTX export."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite

from doublet_rate.io_counts import export_mtx, sanitize_sparse_indices


def _csr_with_indices(data, indices, indptr, shape):
    return sparse.csr_matrix((np.asarray(data, dtype=float), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int32)), shape=shape)


def test_one_based_columns_are_shifted():
    X = _csr_with_indices([10.0, 20.0, 30.0], [1, 2, 3], [0, 2, 3], (2, 3))
    assert int(X.indices.max()) == 3
    Y = sanitize_sparse_indices(X)
    assert Y.shape == (2, 3)
    assert int(Y.indices.max()) == 2
    assert Y[0, 0] == 10.0
    assert Y[0, 1] == 20.0
    assert Y[1, 2] == 30.0


def test_zero_based_oob_entries_are_dropped():
    X = _csr_with_indices([1.0, 9.0], [0, 3], [0, 2, 2], (2, 3))
    Y = sanitize_sparse_indices(X)
    assert Y.shape == (2, 3)
    assert Y.nnz == 1
    assert Y[0, 0] == 1.0


def test_export_mtx_accepts_one_based_csr():
    X = _csr_with_indices([10.0, 20.0, 30.0], [1, 2, 3], [0, 2, 3], (2, 3))

    class _Adata:
        pass

    adata = _Adata()
    adata.X = sanitize_sparse_indices(X)
    adata.obs_names = np.array(["c1", "c2"])
    adata.var_names = np.array(["g1", "g2", "g3"])
    with tempfile.TemporaryDirectory() as td:
        out = export_mtx(adata, Path(td))
        mm = sparse.csc_matrix(mmread(out / "matrix.mtx"))
        assert mm.shape == (3, 2)
        dense = np.asarray(mm.todense())
        assert dense[0, 0] == 10.0
        assert dense[1, 0] == 20.0
        assert dense[2, 1] == 30.0


def test_mmwrite_fails_without_sanitize():
    X = _csr_with_indices([10.0, 20.0, 30.0], [1, 2, 3], [0, 2, 3], (2, 3))
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / "matrix.mtx"
        try:
            mmwrite(path, X.T.tocsc())
        except ValueError as exc:
            assert "exceeds matrix dimension" in str(exc)
            return
        raise AssertionError("expected mmwrite to reject 1-based indices")


if __name__ == "__main__":
    test_one_based_columns_are_shifted()
    test_zero_based_oob_entries_are_dropped()
    test_export_mtx_accepts_one_based_csr()
    test_mmwrite_fails_without_sanitize()
    print("ok")
