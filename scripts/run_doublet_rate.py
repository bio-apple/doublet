#!/usr/bin/env python3
"""Run RNA-only doublet scoring, data-driven Calls, and Predicted Doublet Rates."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from calls import griffiths_mad_calls, rate_from_calls  # noqa: E402
from io_counts import (  # noqa: E402
    SIZE_GATE,
    InputError,
    export_mtx,
    load_sample,
    sample_stem,
)

ALWAYS = ["scdblfinder", "scrublet", "cxds", "bcds", "hybrid"]
GATED = ["doubletdetection", "doubletfinder", "solo"]
NATIVE = {
    "scdblfinder": "native:scdblfinder-dbr.sd=1",
    "scrublet": "native:scrublet-bimodality",
}


def _read_skip(det_dir: Path, name: str) -> str | None:
    p = det_dir / f"{name}.skip.txt"
    if p.exists():
        return p.read_text().strip()
    return None


def _read_detector(det_dir: Path, name: str) -> pd.DataFrame | None:
    p = det_dir / f"{name}.tsv"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    if "barcode" not in df.columns or "score" not in df.columns:
        return None
    df["barcode"] = df["barcode"].astype(str)
    return df


def apply_calls(df: pd.DataFrame, name: str) -> tuple[pd.DataFrame, str]:
    native_col = "native_call" if "native_call" in df.columns else None
    if native_col and name in NATIVE:
        calls = df[native_col].astype(str).where(df[native_col].isin(["singlet", "doublet"]), "")
        if (calls != "").any():
            df = df.copy()
            df["call"] = calls
            return df, NATIVE[name]
    df = df.copy()
    df["call"] = griffiths_mad_calls(df["score"].to_numpy())
    return df, "mad-griffiths"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="10x MTX directory, 10x .h5, or single-sample .h5ad")
    parser.add_argument("--output-dir", default=None, help="directory for TSV outputs (default: beside input)")
    parser.add_argument(
        "--fast",
        action="store_true",
        help="skip DoubletDetection, DoubletFinder, and Solo (smoke test / Quick Start)",
    )
    args = parser.parse_args()

    in_path = Path(args.input).expanduser().resolve()
    try:
        adata = load_sample(in_path)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    stem = sample_stem(in_path)
    out_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else in_path.parent
    if in_path.is_dir() and args.output_dir is None:
        out_dir = in_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    n_input = int(adata.n_obs)
    barcodes = pd.Index(adata.obs_names.astype(str), name="barcode")
    run_gated = n_input <= SIZE_GATE and not args.fast
    gated_skip = "fast:skip_gated" if args.fast else f"size_gate:>{SIZE_GATE}"

    work = out_dir / f".{stem}_doublet_work"
    mtx_dir = export_mtx(adata, work / "mtx")
    det_dir = work / "detectors"
    det_dir.mkdir(parents=True, exist_ok=True)

    from detectors_python import run_doubletdetection, run_scrublet, run_solo

    run_scrublet(adata, det_dir)
    if run_gated:
        run_doubletdetection(adata, det_dir)
        run_solo(adata, det_dir)
    else:
        for name in ("doubletdetection", "solo"):
            (det_dir / f"{name}.skip.txt").write_text(gated_skip + "\n")

    rscript = ROOT / "detectors_r.R"
    r_cmd = ["Rscript", str(rscript), "--mtx-dir", str(mtx_dir), "--outdir", str(det_dir)]
    if run_gated:
        r_cmd.append("--doubletfinder")
    else:
        (det_dir / "doubletfinder.skip.txt").write_text(gated_skip + "\n")
    try:
        r = subprocess.run(r_cmd, check=False, capture_output=True, text=True)
    except FileNotFoundError:
        r = subprocess.CompletedProcess(r_cmd, 1, stdout="", stderr="Rscript not found")
    if r.returncode != 0:
        reason = (r.stderr or r.stdout or "Rscript failed").strip().splitlines()
        tail = reason[-1] if reason else "Rscript failed"
        for name in ["scdblfinder", "cxds", "bcds", "hybrid"] + (["doubletfinder"] if run_gated else []):
            skip = det_dir / f"{name}.skip.txt"
            tsv = det_dir / f"{name}.tsv"
            if not skip.exists() and not tsv.exists():
                skip.write_text(f"error:{tail}\n")

    cell = pd.DataFrame({"barcode": barcodes})
    rows = []
    for name in ALWAYS + GATED:
        skip = _read_skip(det_dir, name)
        det = _read_detector(det_dir, name)
        rec = {
            "detector": name,
            "n_input": n_input,
            "n_scored": 0,
            "n_called": 0,
            "n_doublet": 0,
            "predicted_doublet_rate": float("nan"),
            "status": "skipped",
            "skipped_reason": "",
            "call_rule": "",
        }
        if name in GATED and not run_gated:
            rec["skipped_reason"] = gated_skip
            rows.append(rec)
            continue
        if skip and det is None:
            rec["skipped_reason"] = skip
            rows.append(rec)
            continue
        if det is None:
            rec["skipped_reason"] = skip or "error:no detector output"
            rows.append(rec)
            continue
        det, rule = apply_calls(det, name)
        merged = cell[["barcode"]].merge(det[["barcode", "score", "call"]], on="barcode", how="left")
        cell[f"{name}_score"] = merged["score"]
        cell[f"{name}_call"] = merged["call"].fillna("")
        n_scored = int(merged["score"].notna().sum())
        rec.update(rate_from_calls(cell[f"{name}_call"], n_scored))
        rec["status"] = "ran"
        rec["call_rule"] = rule
        rec["skipped_reason"] = ""
        rows.append(rec)

    sample = pd.DataFrame(rows)
    cells_path = out_dir / f"{stem}.doublet_cells.tsv"
    sample_path = out_dir / f"{stem}.doublet_sample.tsv"
    cell.to_csv(cells_path, sep="\t", index=False)
    sample.to_csv(sample_path, sep="\t", index=False)
    print(cells_path)
    print(sample_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
