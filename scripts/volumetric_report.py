"""Volumetric benchmark report: QED-C-style width x depth tables.

Reads the fresh bench_results.json (QASMBench) and head_to_head /
suppress results, emits a volumetric report showing each circuit's
(width, depth, fidelity) position - the format QED-C uses for
cross-platform comparison.  Output: volumetric_report.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def main() -> int:
    bench = json.load(open(REPO / "bench_results.json", encoding="utf-8"))
    rows = [r for r in bench.get("circuits", [])
            if "compactq" in r and "two_qubit" in r.get("compactq", {})]
    md = ["# Volumetric Benchmark Report", "",
          "Width = qubits, Depth = compactq output depth.", ""]
    md.append("| width | depth | circuit | compactq 2q | qiskit 2q |")
    md.append("|---|---|---|---|---|")
    for r in sorted(rows, key=lambda x: (x.get("num_qubits", 0), x["compactq"]["depth"])):
        nq = r.get("num_qubits", "?")
        md.append(f"| {nq} | {r['compactq']['depth']} | {r['circuit']} "
                  f"| {r['compactq']['two_qubit']} "
                  f"| {r.get('qiskit_l3', {}).get('two_qubit', '-')} |")
    md.append("")
    (REPO / "volumetric_report.md").write_text("\n".join(md) + "\n",
                                               encoding="utf-8")
    print(f"volumetric_report.md: {len(rows)} circuits")
    return 0


if __name__ == "__main__":
    sys.exit(main())
