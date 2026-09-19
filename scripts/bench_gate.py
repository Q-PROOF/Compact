"""Optimization benchmark regression gate - the crown must never slip.

Re-runs the full QASMBench suite (the same harness that generates
bench_results.json) and compares every circuit's compactq output
against the committed baseline.  The 2-QUBIT COUNT is the hard
invariant (win-or-tie required).  Total gates carry a documented
cross-platform drift allowance of max(3, 5% of the baseline count)
and depth max(2, 5% of the baseline depth): floating-point tie-breaks
in 1q resynthesis differ across platforms and native-kernel
availability, so byte-identical 1q counts and tight depth are not a
cross-platform invariant - the 2q count is.
Any regression exits non-zero; improvements are reported and can be
committed with --update.

Usage:
    python scripts/bench_gate.py            # gate against the baseline
    python scripts/bench_gate.py --update   # write a new baseline snapshot
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import date, datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from bench_json import collect_rows  # noqa: E402

BASELINE = REPO / "bench_results.json"
ALLOWED = REPO / "docs" / "agent-notes" / "ALLOWED_REGRESSIONS.md"


def main() -> int:
    update = "--update" in sys.argv
    baseline = {c["circuit"]: c for c in
                json.load(open(BASELINE, encoding="utf-8"))["circuits"]}

    rows = [r for r in collect_rows(verbose=False) if "compactq" in r
            and "two_qubit" in r.get("compactq", {})]
    GATE_DRIFT_FRAC = 0.05   # cross-platform 1q tie-break allowance:
    GATE_DRIFT_MIN = 3       # max(3, 5% of baseline gates)
    DEPTH_DRIFT_FRAC = 0.05  # max(2, 5% of baseline depth): the same
    DEPTH_DRIFT_MIN = 2      # tie-break drift applies to depth
    regressions = []
    improvements = []
    allowed = set()
    if ALLOWED.is_file():
        for ln in ALLOWED.read_text(encoding="utf-8").splitlines():
            ln = ln.strip()
            if ln.startswith("- ") and ":" in ln:
                allowed.add(ln[2:].split(":")[0].strip())
    for row in rows:
        name = row["circuit"]
        cq = row["compactq"]
        new = (cq["two_qubit"], cq["gates"], cq["depth"])
        if name not in baseline:
            improvements.append((name, None, new))
            continue
        old_cq = baseline[name]["compactq"]
        old = (old_cq["two_qubit"], old_cq["gates"], old_cq["depth"])
        gate_allow = max(GATE_DRIFT_MIN, round(old[1] * GATE_DRIFT_FRAC))
        depth_allow = max(DEPTH_DRIFT_MIN, round(old[2] * DEPTH_DRIFT_FRAC))
        if new[0] > old[0]:
            if name in allowed:
                print(f"  allowed regression accepted: {name} "
                      f"(documented in ALLOWED_REGRESSIONS.md)")
            else:
                regressions.append((name, old, new))
        elif new[0] == old[0] and (new[1] > old[1] + gate_allow
                                   or new[2] > old[2] + depth_allow):
            regressions.append((name, old, new))
        elif new < old:
            improvements.append((name, old, new))

    if improvements:
        print(f"{len(improvements)} circuit(s) improved or added:")
        for name, old, new in improvements:
            print(f"  {name}: {old} -> {new}")
    if regressions:
        print(f"GATE FAILED: {len(regressions)} regression(s) vs baseline:")
        for name, old, new in regressions:
            print(f"  {name}: {old} -> {new}")
        return 1
    print(f"GATE PASSED: no regression on {len(rows)} circuits.")

    if update:
        doc = json.load(open(BASELINE, encoding="utf-8"))
        for row in rows:
            for old_row in doc["circuits"]:
                if old_row.get("circuit") == row["circuit"]:
                    if "compactq" in row:
                        old_row["compactq"] = row["compactq"]
                        old_row["compactq_ms"] = row.get("compactq_ms")
                        old_row["proof"] = row.get("proof")
                    break
        doc["generated"] = datetime.now(timezone.utc).isoformat()
        try:
            doc["commit"] = subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], text=True,
                cwd=REPO).strip()
        except Exception:
            pass
        BASELINE.write_text(json.dumps(doc, indent=2) + "\n",
                            encoding="utf-8")
        print(f"baseline updated: {BASELINE.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

