"""Machine-readable benchmark harness.

Runs the QASMBench small suite (Compact vs Qiskit L3) and emits:
  bench_results.json  - full per-circuit records: counts, wall time,
                        referee fidelity, proof status, tool versions,
                        commit SHA, timestamp, machine info
  bench_results.md    - generated Markdown summary (never hand-edited)

Usage:
  python scripts/bench_json.py [--max-qubits N] [--out-dir DIR]

Verification outcome per row: "exact-unitary" (dense referee at <=8q),
"unverified" (>8q), or "INEQUIVALENT" (release blocker in the making).
"""
from __future__ import annotations

import json
import platform
import re
import subprocess
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np  # noqa: E402
from qiskit.quantum_info import Operator  # noqa: E402


def git_sha():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, cwd=str(REPO)).stdout.strip()[:12]
    except Exception:
        return "unknown"


def _native_info():
    try:
        import compactq_native
        import importlib.metadata
        return {"installed": True,
                "version": importlib.metadata.version("compactq-native"),
                "kernels": sorted(k for k in dir(compactq_native)
                                  if not k.startswith("_"))}
    except Exception:
        return {"installed": False}


def _markdown(doc):
    lines = [
        "# Benchmark results (generated - do not edit)",
        "",
        f"Generated {doc['generated']} | commit `{doc['commit']}` | "
        f"compactq {doc['package_version']} | qiskit {doc['qiskit']} | "
        f"native {'yes ' + doc['native_kernels'].get('version', '') if doc['native_kernels'].get('installed') else 'no'}",
        "",
        "| circuit | q | Compact T/2q/depth | Qiskit L3 T/2q/depth | proof | referee fid |",
        "|---|---|---|---|---|---|",
    ]
    for r in doc["circuits"]:
        if r.get("status") == "skipped":
            lines.append(f"| {r['circuit']} | - | skipped: {r.get('reason', '')} "
                         "| - | - | - |")
            continue
        cq, qk = r.get("compactq", {}), r.get("qiskit_l3", {})

        def fmt(d):
            if "error" in d:
                return "ERROR"
            return f"{d['gates']}/{d['two_qubit']}/{d['depth']}"

        lines.append(
            f"| {r['circuit']} | {r['num_qubits']} | {fmt(cq)} | {fmt(qk)} "
            f"| {r.get('proof', '-')} | {r.get('referee_fidelity', '-')} |")
    return "\n".join(lines) + "\n"


def collect_rows(max_q: int = 32, verbose: bool = True):
    """Run the QASMBench small suite; returns the per-circuit rows list
    (the same schema bench_results.json carries).  Reused by
    scripts/bench_gate.py for the optimization regression gate."""
    from realbench import (prepare_unitary, normalize_qiskit,
                           normalize_harder, run_qiskit, metrics)
    root = REPO / "third_party" / "QASMBench" / "small"
    files = []
    for d in sorted(root.iterdir()):
        if d.is_dir():
            for qf in sorted(d.glob("*.qasm")):
                if "_transpiled" not in qf.name and "_out" not in qf.name:
                    files.append(qf)

    rows = []
    for f in files:
        name = f.stem
        text = f.read_text()
        if re.search(r"^\s*if\s*\(", text, re.M):
            rows.append({"circuit": name, "status": "skipped",
                         "reason": "classical control flow"})
            continue
        prepared, reason = prepare_unitary(text)
        if prepared is None:
            rows.append({"circuit": name, "status": "skipped", "reason": reason})
            continue
        m = re.search(r"_n(\d+)", f.name)
        nq = int(m.group(1)) if m else prepared.num_qubits
        if nq > max_q:
            rows.append({"circuit": name, "status": "skipped",
                         "reason": f">{max_q} qubits"})
            continue
        try:
            inp = normalize_qiskit(prepared)
        except Exception:
            inp = normalize_harder(prepared)
        row = {"circuit": name, "num_qubits": nq,
               "input": dict(zip(("gates", "two_qubit", "depth"), metrics(inp)))}
        try:
            from compactq.qiskit_bridge import from_qiskit, to_qiskit
            from compactq import optimize_search
            circ = from_qiskit(inp)
            t0 = time.perf_counter()
            opt = optimize_search(circ)
            dt = time.perf_counter() - t0
            out_qc = to_qiskit(opt)
            row["compactq"] = dict(zip(("gates", "two_qubit", "depth"),
                                       metrics(out_qc)))
            row["compactq_ms"] = round(dt * 1000)
            row["proof"] = "exact-unitary" if nq <= 8 else "unverified"
            if nq <= 8:
                ref = Operator(inp).data
                got = Operator(out_qc).data
                fv = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                row["referee_fidelity"] = round(fv, 12)
                if fv < 1 - 1e-6:
                    row["proof"] = "INEQUIVALENT"
        except Exception as e:
            row["compactq"] = {"error": f"{type(e).__name__}: {e}"}
        try:
            qk, dtk = run_qiskit(inp)
            row["qiskit_l3"] = dict(zip(("gates", "two_qubit", "depth"),
                                        metrics(qk)))
            row["qiskit_ms"] = round(dtk * 1000)
            if nq <= 8:
                ref = Operator(inp).data
                got = Operator(qk).data
                fv = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                row["qiskit_referee_fidelity"] = round(fv, 12)
        except Exception as e:
            row["qiskit_l3"] = {"error": f"{type(e).__name__}: {e}"}
        rows.append(row)
        if verbose:
            cq = row.get("compactq", {})
            qkr = row.get("qiskit_l3", {})
            print(f"{name:22s} {cq.get('gates', '-'):>5}/{cq.get('two_qubit', '-'):>4}"
                  f" | qk {qkr.get('gates', '-'):>5}/{qkr.get('two_qubit', '-'):>4}",
                  flush=True)
    return rows


def main() -> int:
    import compactq
    try:
        import qiskit
        qiskit_version = qiskit.__version__
    except Exception:
        qiskit_version = "not installed"

    max_q = 32
    if "--max-qubits" in sys.argv:
        max_q = int(sys.argv[sys.argv.index("--max-qubits") + 1])
    out_dir = REPO
    if "--out-dir" in sys.argv:
        out_dir = Path(sys.argv[sys.argv.index("--out-dir") + 1])
        out_dir.mkdir(exist_ok=True)

    t_start = time.perf_counter()
    rows = collect_rows(max_q=max_q)

    doc = {
        "schema": "qproof-bench/1",
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "package_version": compactq.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "qiskit": qiskit_version,
        "native_kernels": _native_info(),
        "wall_time_s": round(time.perf_counter() - t_start, 1),
        "circuits": rows,
    }
    (Path(out_dir) / "bench_results.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    (Path(out_dir) / "bench_results.md").write_text(
        _markdown(doc), encoding="utf-8")
    print(f"\nwrote bench_results.json + bench_results.md ({len(rows)} rows, "
          f"{doc['wall_time_s']}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
