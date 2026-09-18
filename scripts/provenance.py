"""Shared provenance block for every benchmark artifact.

Every generator embeds `environment()` so each artifact records the exact
git SHA, package version, dependency versions, Python, OS, CPU and
timestamp it was produced with — the evidence-integrity contract.
"""
from __future__ import annotations

import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def git_sha() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, cwd=str(REPO)).stdout.strip()[:12]
    except Exception:
        return "unknown"


def _version(mod_name: str) -> str:
    try:
        mod = __import__(mod_name)
        return str(getattr(mod, "__version__", "unknown"))
    except Exception:
        return "not installed"


def environment(extra_deps: tuple = ()) -> dict:
    deps = {"qiskit": _version("qiskit"),
            "numpy": _version("numpy")}
    for name in extra_deps:
        deps[name] = _version(name)
    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_sha": git_sha(),
        "python": sys.version.split()[0],
        "os": f"{platform.system()} {platform.release()} {platform.machine()}",
        "cpu": platform.processor() or platform.machine(),
        "dependencies": deps,
    }
