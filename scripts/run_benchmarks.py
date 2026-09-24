"""One-command benchmark regeneration — the v0.2.4 reproduction contract.

A stranger with a clean checkout runs:

    python scripts/run_benchmarks.py                # everything
    python scripts/run_benchmarks.py --only repro   # one artifact family

and gets the committed results/ artifacts regenerated from pinned
inputs: the suite is validated against results/MANIFEST.json before any
number is produced (a referee/tool/suite drift aborts loudly), and every
artifact embeds its own environment block (commit, versions, OS).

Families:
  repro       same-input same-gate-set comparison vs Qiskit L3 /
              pytket / Cirq (+ PyZX / BQSKit when installed)
              -> results/repro.json|csv|md   (needs qiskit)
  ddceiling   DD-prover structured-circuit proof ceiling + before/after
              -> results/dd_ceiling.json|md  (zero dependencies)
  latency     import time + small-circuit optimize latency
              -> results/latency.json        (needs qiskit for the
                                              qiskit comparison rows)

QASMBench / MQT Bench runs keep their dedicated scripts
(scripts/mqtbench_run.py, scripts/realbench.py) — they need network
downloads and are documented in results/README.md.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"

# every argument forwarded to a family's main() is chosen from these
# fixed tables — nothing free-form is passed through
FAMILIES = {
    "repro": ("repro_harness.py", "qiskit"),
    "ddceiling": ("dd_ceiling_bench.py", None),
    "latency": ("latency_bench.py", "qiskit"),
}
KNOWN_TOOLS = {"compact", "qiskit", "pytket", "cirq", "pyzx", "bqskit"}


def _have(mod: str) -> bool:
    try:
        __import__(mod)
        return True
    except Exception:
        return False


def main() -> int:
    ap = argparse.ArgumentParser(description="regenerate results/ artifacts")
    ap.add_argument("--only", default=",".join(FAMILIES),
                    help="comma-separated subset: " + ", ".join(FAMILIES))
    ap.add_argument("--tools", default="compact,qiskit,pytket,cirq",
                    help="tools for the repro family")
    args = ap.parse_args()
    wanted = [t.strip() for t in args.only.split(",") if t.strip()]
    unknown = [t for t in wanted if t not in FAMILIES]
    if unknown:
        print(f"unknown families {unknown}; expected {sorted(FAMILIES)}")
        return 1
    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    if not tools or any(t not in KNOWN_TOOLS for t in tools):
        print(f"invalid --tools; expected a comma list from "
              f"{sorted(KNOWN_TOOLS)}")
        return 1

    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))

    rc = 0
    for name in wanted:
        script, dep = FAMILIES[name]
        if dep and not _have(dep):
            print(f"== {name}: SKIPPED ({dep} not installed; "
                  f"pip install '{dep}') ==")
            continue
        mod = importlib.import_module(script[:-3])
        argv = [script]
        if name == "repro":
            argv += ["--tools", ",".join(tools)]
        print(f"== {name}: {script} {' '.join(argv[1:])} ==", flush=True)
        saved = sys.argv
        try:
            sys.argv = argv
            ret = mod.main()
        finally:
            sys.argv = saved
        if ret:
            print(f"== {name}: FAILED (exit {ret}) ==")
            rc = ret or 1
        else:
            print(f"== {name}: OK ==")
    if rc:
        print("\none or more benchmark families FAILED — results/ is not "
              "consistent; do not publish.")
    else:
        print("\nall requested benchmark families regenerated.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
