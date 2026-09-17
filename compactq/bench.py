"""Benchmark CLI: compactq vs. raw circuits (and vs. Qiskit when installed).

Usage:
    python -m compactq.bench [--quick] [--out FILE]

Fairness rules: every optimizer receives the *same* input circuit (our QASM
serialization), the same basis set, fixed seeds, and best-of-3 timings.
"""
from __future__ import annotations

import argparse
import time

from .benchmarks import default_suite
from .io_qasm import from_qasm, to_qasm
from .optimize import optimize


def _two_q(ops) -> int:
    return sum(1 for g in ops if len(g.qubits) == 2)


def _metrics(circ):
    return len(circ.ops), _two_q(circ.ops), circ.depth()


def _qiskit_baseline(qasm: str):
    """Run Qiskit optimization_level=3 on the same QASM. Returns (gates, 2q, depth, ms) or None."""
    try:
        from qiskit import transpile
        from qiskit import qasm2 as qk_qasm
    except Exception:
        return None
    qc = qk_qasm.loads(qasm)
    basis = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p", "cx", "cz", "swap"]
    best = None
    for _ in range(3):
        t0 = time.perf_counter()
        out = transpile(qc, basis_gates=basis, optimization_level=3, seed_transpiler=42)
        dt = (time.perf_counter() - t0) * 1000.0
        best = dt if best is None else min(best, dt)
    ops = out.count_ops()
    total = int(sum(v for k, v in ops.items() if k not in ("barrier",)))
    two_q = int(ops.get("cx", 0) + ops.get("cz", 0) + ops.get("swap", 0)
                + ops.get("ecr", 0) + ops.get("rxx", 0))
    depth = int(out.depth())
    return total, two_q, depth, best


def run_suite(quick: bool = False):
    suite = default_suite()
    if quick:
        suite = {k: v for k, v in suite.items() if k in ("ghz-5", "qft-3", "clifford-ladder-4")}
    rows = []
    for name, circ in suite.items():
        qasm = to_qasm(circ)
        raw = _metrics(circ)

        best_ms = None
        for _ in range(3):
            t0 = time.perf_counter()
            opt = optimize(from_qasm(qasm), verify=circ.num_qubits <= 6)
            dt = (time.perf_counter() - t0) * 1000.0
            best_ms = dt if best_ms is None else min(best_ms, dt)
        ours = _metrics(opt)

        qk = _qiskit_baseline(qasm)
        rows.append((name, raw, ours, best_ms, qk))
    return rows


def _fmt_row(name, raw, ours, ours_ms, qk):
    def cell(v):
        return str(v) if v is not None else "n/a"
    line = f"| {name} | {raw[0]} / {raw[1]} / {raw[2]} | {ours[0]} / {ours[1]} / {ours[2]} | {ours_ms:.2f} |"
    if qk is not None:
        qt, q2, qd, qms = qk
        gain = f"-{(1 - ours[1] / q2) * 100:.0f}%" if q2 else "n/a"
        line += f" {qt} / {q2} / {qd} | {qms:.2f} | {gain} |"
    return line


def main() -> int:
    ap = argparse.ArgumentParser(prog="compactq-bench", description=__doc__)
    ap.add_argument("--quick", action="store_true", help="run a 3-circuit quick suite")
    ap.add_argument("--out", default=None, help="also write the markdown table to FILE")
    args = ap.parse_args()

    rows = run_suite(quick=args.quick)
    has_qk = rows and rows[0][4] is not None

    lines = []
    header = "| circuit | raw gates/2q/depth | compactq gates/2q/depth | compactq ms |"
    sep = "|---|---|---|---|"
    if has_qk:
        header += " qiskit L3 gates/2q/depth | qiskit ms | 2q-gate gain vs qiskit |"
        sep += "---|---|---|"
    lines.append(header)
    lines.append(sep)
    for name, raw, ours, ms, qk in rows:
        lines.append(_fmt_row(name, raw, ours, ms, qk))

    table = "\n".join(lines)
    print(table)
    if args.out:
        from pathlib import Path
        out_path = Path(args.out).expanduser().resolve()
        if not out_path.parent.is_dir():
            raise SystemExit(f"--out directory does not exist: "
                             f"{out_path.parent}")
        out_path.write_text(table + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
