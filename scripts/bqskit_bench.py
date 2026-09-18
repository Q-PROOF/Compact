"""BQSKit column for the benchmark matrix (benchmark suite v1).

Runs BQSKit's documented compile pipeline (`bqskit.compile`) on the same
normalized QASMBench small circuits as scripts/bench_json.py, referees
every output with qiskit's Operator, and writes results/bqskit.json in
the standard artifact schema.

    pip install bqskit
    python scripts/bqskit_bench.py

Honesty rule (contributor rule #2): no bqskit numbers may appear in the
README until produced by this script.
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
warnings.filterwarnings("ignore")

try:
    import bqskit
    from bqskit import compile as bqskit_compile
except ImportError:
    print("bqskit is not installed - nothing to measure.")
    print("  pip install bqskit && python scripts/bqskit_bench.py")
    sys.exit(0)


def _cx_eq(qc):
    """Level-B metric for a qiskit circuit (SWAP = 3 CX, others = 1)."""
    return sum(3 if i.operation.name == "swap" else 1
               for i in qc.data if len(i.qubits) == 2)


def main() -> int:
    import numpy as np  # noqa: F401
    from qiskit import qasm2
    from qiskit.quantum_info import Operator

    from bench_json import git_sha
    from realbench import metrics, normalize_qiskit, prepare_unitary

    records = []
    root = REPO / "third_party" / "QASMBench" / "small"
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        for qf in sorted(d.glob("*.qasm")):
            if "_transpiled" in qf.name or "_out" in qf.name:
                continue
            prepared, _reason = prepare_unitary(qf.read_text())
            if prepared is None:
                continue
            try:
                inp = normalize_qiskit(prepared)
            except Exception:
                continue
            ref = Operator(inp).data
            try:
                t0 = time.perf_counter()
                tmp_qasm = REPO / "results" / f"_bqskit_{qf.stem}.qasm"
                tmp_qasm.write_text(qasm2.dumps(inp), encoding="utf-8")
                bqs = bqskit.Circuit.from_file(str(tmp_qasm))
                compiled = bqskit_compile(bqs)
                out_path = REPO / "results" / f"_bqskit_out_{qf.stem}.qasm"
                compiled.save(str(out_path))
                dt = time.perf_counter() - t0
                out_qc = qasm2.loads(out_path.read_text(encoding="utf-8"))
                out_path.unlink(missing_ok=True)
                tmp_qasm.unlink(missing_ok=True)
                got = Operator(out_qc).data
                fid = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                records.append({
                    "circuit": qf.stem,
                    "tool": "bqskit",
                    "bqskit_version": getattr(bqskit, "__version__",
                                              "unknown"),
                    "num_qubits": inp.num_qubits,
                    "gates": sum(1 for _ in out_qc.data),
                    "two_qubit": sum(1 for i in out_qc.data
                                     if len(i.qubits) == 2),
                    "cx_equivalent": _cx_eq(out_qc),
                    "depth": out_qc.depth(),
                    "output_fidelity": round(fid, 12),
                    "status": "OK" if fid > 1 - 1e-6 else "INEQUIVALENT",
                    "wall_ms": round(dt * 1000),
                })
                print(f"{qf.stem:24s} gates {records[-1]['gates']:5d} "
                      f"2q {records[-1]['two_qubit']:4d} "
                      f"{records[-1]['status']}")
            except Exception as e:
                records.append({"circuit": qf.stem, "tool": "bqskit",
                                "status": f"error {type(e).__name__}: {e}"})
                print(f"{qf.stem:24s} ERROR {type(e).__name__}: {e}")

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "referee": "qiskit.quantum_info.Operator, |Tr(U+V)|/d > 1 - 1e-6",
        "pipeline": "bqskit.compile (documented default optimization)",
        "count": len(records),
        "records": records,
    }
    out = REPO / "results" / "bqskit.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(records)} records)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
