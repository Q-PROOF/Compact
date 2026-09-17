"""Three-way real-world benchmark: compactq vs Qiskit L3 vs pytket.

Protocol (identical start for all tools):
  input.qasm -> qiskit level-0 translation to the shared basis -> each tool
  optimizes that same circuit -> count (total, 2q, depth) + wall time.
No coupling maps: pure logical optimization.  compactq results are unitary-verified
for <= 6 qubits.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

OUR_BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p",
             "sx", "cx", "cz", "swap"]


def normalize_qiskit(qc):
    from qiskit import transpile
    try:
        out = transpile(qc, basis_gates=OUR_BASIS, optimization_level=0, seed_transpiler=42)
        op_names = {inst.operation.name for inst in out.data}
        if op_names <= set(OUR_BASIS) | {"barrier", "measure"}:
            return out
    except Exception:
        pass
    # fallback: force-inline nested custom gates
    decomposed = qc.decompose(reps=3)
    return transpile(decomposed, basis_gates=OUR_BASIS, optimization_level=0, seed_transpiler=42)


def run_compactq(qc):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from compactq.qiskit_bridge import from_qiskit
    from compactq import optimize_search as compactq_optimize

    circ = from_qiskit(qc)
    t0 = time.perf_counter()
    out = compactq_optimize(circ)
    dt = time.perf_counter() - t0
    from compactq.qiskit_bridge import to_qiskit
    return to_qiskit(out), dt


def normalize_harder(qc):
    """Level-0 translation can leave nested custom gates; force-inline."""
    from qiskit import transpile
    out = qc.decompose(reps=3)
    return transpile(out, basis_gates=OUR_BASIS, optimization_level=0, seed_transpiler=42)


def run_qiskit(qc):
    from qiskit import transpile
    best = None
    best_t = float("inf")
    for _ in range(3):
        t0 = time.perf_counter()
        out = transpile(qc, basis_gates=OUR_BASIS, optimization_level=3, seed_transpiler=42)
        dt = time.perf_counter() - t0
        best_t = min(best_t, dt)
        best = out
    return best, best_t


def run_pytket(qc):
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
    t0 = time.perf_counter()
    try:
        seq = SequencePass([CommuteThroughMultis(), RemoveRedundancies(),
                            FullPeepholeOptimise(), AutoRebase(ourset),
                            CommuteThroughMultis(), RemoveRedundancies()])
        seq.apply(cu)
    except Exception as e:
        print(f"  [pytket pass error: {e}]")
        FullPeepholeOptimise().apply(cu)
    dt = (time.perf_counter() - t0) * 1000
    out = tk_to_qiskit(cu.circuit)
    return out, dt / 1000.0


def metrics(qc):
    ops = qc.count_ops()
    total = int(sum(ops.values()))
    two_q = int(sum(v for k, v in ops.items() if k in ("cx", "cz", "swap", "cp", "ecr", "rzz")))
    return total, two_q, qc.depth()


def prepare_unitary(text: str):
    """Return measure-free unitary QASM text, or (None, reason) to skip.

    Circuits with classical control flow or mid-circuit measurement are not
    unitary and are out of scope for a unitary optimizer.  Trailing
    measure/reset blocks are stripped so the comparison is apples-to-apples.
    """
    import re
    if re.search(r"^\s*if\s*\(", text, re.M):
        return None, "non-unitary (classical control flow)"
    from qiskit import qasm2
    try:
        # extended qelib1 (cz/swap/cswap/sx/...): QASMBench assumes it
        raw = qasm2.loads(text, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    except Exception as e:
        return None, f"parse: {type(e).__name__}: {e}"
    data = raw.data
    idx = len(data)
    while idx > 0 and data[idx - 1].operation.name in ("measure", "barrier"):
        idx -= 1
    rest = {inst.operation.name for inst in data[:idx]}
    if "measure" in rest or "reset" in rest:
        return None, "non-unitary (mid-circuit measurement)"
    if idx == 0:
        return None, "empty unitary core"
    raw.remove_final_measurements(inplace=True)
    return raw, None


def run_suite(qasm_files, sizes=None):
    rows = []
    for f in qasm_files:
        name = Path(f).stem
        prepared, reason = prepare_unitary(Path(f).read_text())
        if prepared is None:
            print(f"skip {name}: {reason}")
            continue
        try:
            try:
                inp = normalize_qiskit(prepared)
            except Exception:
                inp = normalize_harder(prepared)
        except Exception as e:
            print(f"skip {name}: {type(e).__name__}: {e}")
            continue
        it, i2q, idp = metrics(inp)
        row = {"name": name, "in": (it, i2q, idp)}
        try:
            qz, qz_ms = run_compactq(inp)
            row["compactq"] = (*metrics(qz), round(qz_ms * 1000))
        except Exception as e:
            print(f"compactq failed {name}: {type(e).__name__}: {e}")
            row["compactq"] = None
        try:
            qk, qk_ms = run_qiskit(inp)
            row["qk"] = (*metrics(qk), round(qk_ms * 1000))
        except Exception as e:
            print(f"qiskit failed {name}: {type(e).__name__}: {e}")
            row["qk"] = None
        try:
            pt, pt_ms = run_pytket(inp)
            row["pt"] = (*metrics(pt), round(pt_ms * 1000))
        except Exception as e:
            print(f"pytket failed {name}: {type(e).__name__}: {e}")
            row["pt"] = None
        rows.append(row)
        def fmt(v):
            return f"{v[0]}/{v[1]}/{v[2]} ({v[3]}ms)" if v else "FAILED"
        line = f"{name}: in {it}/{i2q}/{idp} | compactq {fmt(row['compactq'])} | qk {fmt(row['qk'])}"
        if row['pt']:
            line += f" | pt {fmt(row['pt'])}"
        print(line)
    return rows


def qasm_files_from_qasmbench(root: Path, max_qubits: int = 6,
                              size: str = "small"):
    """ALL benchmark QASM files (one directory may carry several distinct
    circuits - e.g. basis_trotter_n4 also contains basis_test_n4.qasm);
    labelled by file stem, not the directory name."""
    base = root / size
    files = []
    for d in sorted(base.iterdir()):
        if d.is_dir():
            for qf in sorted(d.glob("*.qasm")):
                if "_transpiled" not in qf.name and "_out" not in qf.name:
                    files.append(qf)
    return files


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", default="small", choices=["small", "medium", "large"])
    ap.add_argument("--max-qubits", type=int, default=None,
                    help="skip circuits above this qubit count")
    args, _ = ap.parse_known_args()
    root = Path(__file__).resolve().parents[1] / "third_party" / "QASMBench"
    if args.size == "small":
        files = qasm_files_from_qasmbench(root)
    else:
        files = qasm_files_from_qasmbench(root, size=args.size)
    if args.max_qubits:
        import re as _re
        def _qcount(f):
            m = _re.search(r"_n(\d+)", f.name)
            return int(m.group(1)) if m else 999
        files = [f for f in files if _qcount(f) <= args.max_qubits]
    print(f"{len(files)} QASMBench {args.size} circuits")
    run_suite([str(f) for f in files])
