"""BQSKit column for the benchmark matrix (benchmark suite v1, roadmap).

Runs BQSKit's synthesis pipeline on the same normalized QASMBench small
circuits as scripts/bench_json.py, referees every output with qiskit's
Operator, and writes results/bqskit.json in the standard artifact schema.

Status: SCAFFOLD.  BQSKit is not yet measured for the README tables; this
script exists so the column lands in one command once BQSKit is installed:

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
    import bqskit  # noqa: F401
except ImportError:
    print("bqskit is not installed - nothing to measure.")
    print("  pip install bqskit && python scripts/bqskit_bench.py")
    sys.exit(0)


def main() -> int:
    import bqskit
    from bqskit.compiler import Compiler
    from bqskit.passes import UnrollPass, VariableLayerization  # noqa: F401
    from qiskit.quantum_info import Operator

    from realbench import prepare_unitary, normalize_qiskit
    from bench_json import git_sha, _cx_eq

    records = []
    root = REPO / "third_party" / "QASMBench" / "small"
    with Compiler() as compiler:
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
                    bqs = bqskit.Circuit.from_qasm(inp.qasm())
                    bqs = compiler.compile(bqs)
                    dt = time.perf_counter() - t0
                    out_qc = bqs.to_qiskit()
                    got = Operator(out_qc).data
                    fid = float(abs((ref.conj().T @ got).trace()) / ref.shape[0])
                    records.append({
                        "circuit": qf.stem,
                        "tool": "bqskit",
                        "bqskit_version": getattr(bqskit, "__version__", "unknown"),
                        "num_qubits": inp.num_qubits,
                        "gates": sum(1 for _ in out_qc.data),
                        "two_qubit": sum(1 for i in out_qc.data if len(i.qubits) == 2),
                        "cx_equivalent": _cx_eq(out_qc),
                        "output_fidelity": round(fid, 12),
                        "status": "OK" if fid > 1 - 1e-6 else "INEQUIVALENT",
                        "wall_ms": round(dt * 1000),
                    })
                except Exception as e:
                    records.append({"circuit": qf.stem, "tool": "bqskit",
                                    "status": f"error {type(e).__name__}: {e}"})

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "referee": "qiskit.quantum_info.Operator, |Tr(U+V)|/d > 1 - 1e-6",
        "count": len(records),
        "records": records,
    }
    out = REPO / "results" / "bqskit.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(records)} records)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
