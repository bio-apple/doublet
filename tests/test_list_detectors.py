"""CLI --list-detectors prints the static roster."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from doublet_rate.run_doublet_rate import ALWAYS, GATED, format_detector_list


def test_format_detector_list():
    text = format_detector_list()
    lines = text.strip().splitlines()
    assert len(lines) == 8
    assert "available" not in text
    assert "DoubletDecon" not in text
    names = [ln.split()[0] for ln in lines]
    roles = [ln.split()[-1] for ln in lines]
    assert names == [
        "scDblFinder",
        "Scrublet",
        "cxds",
        "bcds",
        "hybrid",
        "DoubletDetection",
        "DoubletFinder",
        "Solo",
    ]
    assert roles == ["always"] * len(ALWAYS) + ["gated"] * len(GATED)


def test_cli_list_detectors_no_input():
    r = subprocess.run(
        [sys.executable, "-m", "doublet_rate", "--list-detectors"],
        cwd=_SRC,
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, r.stderr
    assert r.stdout == format_detector_list()


if __name__ == "__main__":
    test_format_detector_list()
    test_cli_list_detectors_no_input()
    print("ok")
