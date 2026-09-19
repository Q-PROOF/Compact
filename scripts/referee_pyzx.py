"""Independent referee #2: PyZX (ZX-calculus) cross-check.

Re-answers the Compact-vs-input equivalence question through a COMPLETELY
independent implementation: both circuits are loaded into PyZX, reduced
to ZX-graph normal form, and compared as tensors.  Agreement between this
referee and the Qiskit-Operator referee (bench_results.json, gauntlet)
is evidence that the verdicts are implementation-independent.

Emits results/pyzx_referee.json.  Errors are recorded as errors — no
claims are made from failed PyZX runs.
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
warnings.filterwarnings("ignore")


def _zx_tensor(qc):
    import pyzx
    from qiskit import qasm2
    buf = io.StringIO()
    qasm2.dump(qc, buf)
    circ = pyzx.Circuit.from_qasm(buf.getvalue())
    graph = circ.to_graph()
    pyzx.full_reduce(graph, quiet=True)
    return graph.to_tensor()


import io  # noqa: E402


def main() -> int:
    import numpy as np
    import pyzx

    import compactq
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from provenance import environment, git_sha
    from realbench import normalize_qiskit, prepare_unitary

    print(f"pyzx {pyzx.__version__} | compactq {compactq.__version__}")
    root = REPO / "third_party" / "QASMBench" / "small"
    records = []
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        qf = d / f"{d.name}.qasm"
        if not qf.is_file():
            continue
        rec = {"circuit": d.name, "tool": "compactq",
               "referee": "pyzx-zx-tensor"}
        prepared, _reason = prepare_unitary(qf.read_text())
        if prepared is None or prepared.num_qubits > 8:
            rec["status"] = "skipped"
            records.append(rec)
            continue
        try:
            inp = normalize_qiskit(prepared)
            circ = from_qiskit(inp)
            opt = compactq.optimize_search(circ)
            t_in = _zx_tensor(inp)
            t_out = _zx_tensor(to_qiskit(opt))
            d2 = t_in.size ** 0.5
            fid = float(abs(np.sum(np.conj(t_in) * t_out)) / d2)
            rec["output_fidelity"] = round(fid, 12)
            rec["status"] = "OK" if fid > 1 - 1e-6 else "INEQUIVALENT"
        except Exception as e:
            rec["status"] = f"error {type(e).__name__}: {e}"
        records.append(rec)
        print(f"{d.name:24s} {rec['status']}", flush=True)

    ok = sum(1 for r in records if r.get("status") == "OK")
    bad = sum(1 for r in records if r.get("status") == "INEQUIVALENT")
    skipped = sum(1 for r in records if r.get("status") == "skipped")
    errs = len(records) - ok - bad - skipped
    print(f"pyzx referee: {ok} OK, {bad} INEQUIVALENT, {skipped} skipped, "
          f"{errs} errors")

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "referee": ("pyzx ZX-graph full-reduction + tensor comparison "
                    "(independent of qiskit's Operator)"),
        "environment": environment(extra_deps=("pyzx",)),
        "count": len(records),
        "ok": ok,
        "inequivalent": bad,
        "records": records,
    }
    out = REPO / "results" / "pyzx_referee.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
