#!/usr/bin/env python3
"""One-click smoke test: load test/ and run the fast Detectors."""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "test"
OUT = ROOT / "test_out"
PY_PKGS = ("scanpy", "numpy", "pandas", "scipy", "scrublet")
R_PKGS = ("SingleCellExperiment", "scDblFinder", "scds")


def _ok(label: str, ok: bool, detail: str = "") -> None:
    mark = "ok" if ok else "FAIL"
    extra = f"  {detail}" if detail else ""
    print(f"[{mark}] {label}{extra}", flush=True)


def check_python() -> bool:
    good = True
    for name in PY_PKGS:
        try:
            importlib.import_module(name)
            _ok(name, True)
        except ImportError as exc:
            _ok(name, False, str(exc))
            good = False
    return good


def check_r() -> bool:
    rscript = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)
    if rscript.returncode != 0:
        _ok("Rscript", False, "not found")
        return False
    _ok("Rscript", True)
    good = True
    for pkg in R_PKGS:
        cmd = ["Rscript", "-e", f"quit(status=if (requireNamespace('{pkg}', quietly=TRUE)) 0 else 1)"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        _ok(pkg, r.returncode == 0, "missing" if r.returncode else "")
        if r.returncode != 0:
            good = False
    return good


def main() -> int:
    print("== environment ==", flush=True)
    py_ok = check_python()
    r_ok = check_r()
    if not py_ok:
        print("install Python deps: python3 -m pip install -r scripts/requirements.txt", file=sys.stderr)
        return 1
    if not DATA.is_dir() or not (DATA / "matrix.mtx").exists():
        print("error: missing test/ 10x MTX (matrix.mtx + barcodes.tsv + genes.tsv)", file=sys.stderr)
        return 1

    print("== run --fast on test/ ==", flush=True)
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_doublet_rate.py"),
        str(DATA),
        "--output-dir",
        str(OUT),
        "--fast",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stdout:
        print(proc.stdout, end="")
    if proc.stderr:
        print(proc.stderr, end="", file=sys.stderr)
    if proc.returncode != 0:
        return proc.returncode

    sample = OUT / "test.doublet_sample.tsv"
    if not sample.exists():
        print("error: missing", sample, file=sys.stderr)
        return 1
    text = sample.read_text()
    print(text, end="")
    n_ran = sum(1 for line in text.splitlines()[1:] if "\tran\t" in line)
    if n_ran == 0:
        print("error: no Detector ran", file=sys.stderr)
        return 1
    if not r_ok:
        print("R packages missing: some Detectors skipped. Rscript scripts/install_r_packages.R")
    print("demo ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
