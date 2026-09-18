"""Referee pytket's outputs on the QASMBench small suite - both modes.

Mirrors scripts/bench_json.py's proven pipeline (prepare_unitary ->
normalize_qiskit -> tool -> Operator referee) with pytket as the tool,
in default AND safe (allow_swaps=False) modes, so the README's pytket
notes cite fresh, same-protocol evidence.
"""
import re
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

warnings.filterwarnings("ignore")


def run_pytket_safe(qc):
    """FullPeepholeOptimise(allow_swaps=False) - pytket's safe mode
    (2.18.1 note: safe mode is itself inequivalent on 2 of 30 QASMBench
    small circuits - basis_trotter_n4 and sat_n7)."""
    import pytket
    from pytket.extensions.qiskit import qiskit_to_tk, tk_to_qiskit
    from pytket.passes import (AutoRebase, CommuteThroughMultis,
                               FullPeepholeOptimise, RemoveRedundancies,
                               SequencePass)
    from pytket.predicates import CompilationUnit
    OT = pytket.circuit.OpType
    tk = qiskit_to_tk(qc)
    ourset = {OT.H, OT.X, OT.Y, OT.Z, OT.S, OT.Sdg, OT.T, OT.Tdg,
              OT.Rx, OT.Ry, OT.Rz, OT.U1, OT.SX, OT.CX, OT.CZ, OT.SWAP}
    cu = CompilationUnit(tk)
    seq = SequencePass([CommuteThroughMultis(), RemoveRedundancies(),
                        FullPeepholeOptimise(allow_swaps=False),
                        AutoRebase(ourset),
                        CommuteThroughMultis(), RemoveRedundancies()])
    seq.apply(cu)
    return tk_to_qiskit(cu.circuit)


def main() -> int:
    import pytket
    from qiskit.quantum_info import Operator
    from realbench import (prepare_unitary, normalize_qiskit, run_pytket,
                           metrics)

    print(f"pytket {pytket.__version__}")
    root = REPO / "third_party" / "QASMBench" / "small"
    rows = []   # (name, mode, nq, fid, verdict, metrics)
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        for qf in sorted(d.glob("*.qasm")):
            if "_transpiled" in qf.name or "_out" in qf.name:
                continue
            name = qf.stem
            text = qf.read_text()
            if re.search(r"^\s*if\s*\(", text, re.M):
                continue
            prepared, reason = prepare_unitary(text)
            if prepared is None:
                continue
            m = re.search(r"_n(\d+)", name)
            nq = int(m.group(1)) if m else prepared.num_qubits
            if nq > 8:
                continue
            try:
                inp = normalize_qiskit(prepared)
            except Exception:
                continue
            ref = Operator(inp).data
            for mode, runner in (("default", "realbench"), ("safe", None)):
                try:
                    if runner == "realbench":
                        out, _dt = run_pytket(inp)
                    else:
                        out = run_pytket_safe(inp)
                    got = Operator(out).data
                    fv = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                    mm = metrics(out)
                except Exception as e:
                    rows.append((name, mode, nq, None,
                                 f"error {type(e).__name__}", None))
                    continue
                rows.append((name, mode, nq, fv,
                             "OK" if fv > 1 - 1e-6 else "INEQUIVALENT", mm))
    for mode in ("default", "safe"):
        sub = [r for r in rows if r[1] == mode]
        bad = [r for r in sub if r[4] != "OK"]
        for name, _mode, nq, fv, note, _met in sub:
            if note != "OK":
                fstr = f"{fv:.6f}" if fv is not None else "n/a"
                print(f"{mode:8s} {name:24s} nq={nq} fid={fstr}  <-- {note}")
        print(f"{mode} mode: {len(bad)} inequivalent of {len(sub)} refereed")

    # committed raw artifact: one JSON record per refereed circuit, so the
    # inequivalence claims are directly reproducible and auditable
    import json
    from datetime import datetime, timezone
    from bench_json import git_sha
    records = []
    for name, mode, nq, fv, verdict, met in rows:
        records.append({
            "circuit": name,
            "tool": "pytket",
            "pytket_version": pytket.__version__,
            "mode": mode,
            "num_qubits": nq,
            "output_fidelity": fv,
            "status": verdict,
            "metrics": met,
        })
    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "referee": "qiskit.quantum_info.Operator, |Tr(U+V)|/d > 1 - 1e-6",
        "environment": __import__("provenance").environment(
            extra_deps=("pytket",)),
        "count": len(records),
        "records": records,
    }
    out = REPO / "results" / "pytket_referee.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out} ({len(records)} records)")

    # 2q tally vs compactq (from the fresh bench_results.json), invalid
    # pytket outputs excluded from the comparison
    import json
    base = {r["circuit"]: r["compactq"]
            for r in json.load(open(REPO / "bench_results.json",
                                    encoding="utf-8"))["circuits"]
            if "compactq" in r}
    for mode in ("default", "safe"):
        w = t = l = 0
        for name, m, _nq, fv, verdict, met in rows:
            if m != mode or verdict != "OK" or name not in base:
                continue
            if met[1] > base[name]["two_qubit"]:
                w += 1
            elif met[1] == base[name]["two_qubit"]:
                t += 1
            else:
                l += 1
        print(f"compactq vs pytket-{mode} (valid outputs only): "
              f"2q win {w} / tie {t} / lose {l}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
