#!/usr/bin/env python3
"""Shim: run the packaged CLI. Prefer `doublet-rate` or `python -m doublet_rate`."""

from __future__ import annotations

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from doublet_rate.run_doublet_rate import main

if __name__ == "__main__":
    raise SystemExit(main())
