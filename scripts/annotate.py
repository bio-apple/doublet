"""Shim for repo checkouts. Prefer: from doublet_rate import detect_doublets."""

from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from doublet_rate.annotate import annotate_anndata, attach_obs, detect_doublets  # noqa: F401
