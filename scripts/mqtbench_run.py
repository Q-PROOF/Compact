"""MQT Bench four-way benchmark: compactq vs Qiskit L3 vs pytket vs Cirq.

Protocol (identical to scripts/realbench.py, extended to more tools):
  MQT Bench algorithm-level circuit -> strip trailing measure/barrier
  (reject mid-circuit measurement) -> level-0 normalization to the shared
  basis -> every tool optimizes that SAME circuit -> count (total, 2q,
  depth) + wall time.  For circuits <= 8 qubits every tool's output is
  refereed with qiskit.quantum_info.Operator (|Tr(A^dag B)|/d); Qiskit's
  output has elided permutations reified first (fairness), pytket's raw
  default-mode output is refereed as-is.

Tools:
  compactq   optimize_search (equivalence proof included in its time)
  qk     qiskit transpile optimization_level=3, best of 3
  pt     pytket AutoRebase + FullPeepholeOptimise pipeline (realbench's)
  cirq   cirq.optimize_for_target_gateset(CZTargetGateset)

staq is not run: it has no PyPI distribution (the PyPI 'staq' package is
an unrelated C-decompiler) and building github.com/softwareqinc/staq
needs a C++17 toolchain this environment does not have.

Extras:  pip install mqt-bench cirq ply
Exit code 0 unless the infrastructure itself fails.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

OUR_BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p",
             "sx", "cx", "cz", "swap"]
REF_QUBITS = 8

# (benchmark, circuit_size) - sizes within the families MQT documents;
# invalid combos are skipped with their exception at runtime.
CANDIDATES = [
    ("ghz", 4), ("ghz", 6), ("ghz", 8),
    ("wstate", 4), ("wstate", 6),
    ("graphstate", 4), ("graphstate", 6), ("graphstate", 8),
    ("qaoa", 4), ("qaoa", 6), ("qaoa", 8),
    ("qft", 4), ("qft", 6), ("qft", 8),
    ("qftentangled", 4), ("qftentangled", 6),
    ("vqe_real_amp", 4), ("vqe_real_amp", 6), ("vqe_real_amp", 8),
    ("vqe_two_local", 4), ("vqe_two_local", 6),
    ("qnn", 4), ("qnn", 6),
    ("qpeexact", 4), ("qpeexact", 6),
    ("qpeinexact", 4), ("qpeinexact", 6),
    ("qwalk", 4), ("qwalk", 6),
    ("randomcircuit", 4), ("randomcircuit", 6), ("randomcircuit", 8),
    ("cdkm_ripple_carry_adder", 4), ("cdkm_ripple_carry_adder", 6),
    ("draper_qft_adder", 4), ("draper_qft_adder", 6),
    ("grover", 4), ("grover", 6),
    ("bv", 5), ("bv", 7),
    ("dj", 5), ("dj", 7),
    ("ae", 4), ("ae", 6),
]


def log(msg):
    print(msg, flush=True)


def gen_unitary_mqt(name: str, size: int):
    """MQT circuit -> (qiskit_circuit, None) or (None, reason)."""
    from mqt.bench import BenchmarkLevel, get_benchmark
    try:
        qc = get_benchmark(name, level=BenchmarkLevel.ALG, circuit_size=size)
    except Exception as e:  # family/size not supported
        return None, f"gen: {type(e).__name__}: {str(e)[:60]}"
    data = qc.data
    idx = len(data)
    while idx > 0 and data[idx - 1].operation.name in ("measure", "barrier"):
        idx -= 1
    rest = {inst.operation.name for inst in data[:idx]}
    if "measure" in rest or "reset" in rest:
        return None, "non-unitary (mid-circuit measurement)"
    if idx == 0:
        return None, "empty unitary core"
    qc.remove_final_measurements(inplace=True)
    return qc, None


def normalize(qc):
    from qiskit import transpile
    out = transpile(qc, basis_gates=OUR_BASIS, optimization_level=0,
                    seed_transpiler=42)
    names = {inst.operation.name for inst in out.data}
    if not names <= set(OUR_BASIS) | {"barrier"}:
        out = transpile(qc.decompose(reps=3), basis_gates=OUR_BASIS,
                        optimization_level=0, seed_transpiler=42)
    return out


def metrics(qc):
    ops = qc.count_ops()
    total = int(sum(ops.values()))
    two_q = int(sum(v for k, v in ops.items() if k in ("cx", "cz", "swap", "cp", "ecr", "rzz")))
    return total, two_q, qc.depth()


def fid(u, v):
    return float(abs(np.sum(np.conj(u) * v)) / u.shape[0])


def reify_perm(qc_out):
    """Reify Qiskit-elided permutations as explicit SWAPs (fair referee)."""
    from qiskit import QuantumCircuit
    if qc_out.layout is None:
        return qc_out, 0
    nq = qc_out.num_qubits
    perm = qc_out.layout.final_index_layout()
    if perm == list(range(nq)):
        return qc_out, 0
    fixed = QuantumCircuit(nq)
    for inst in qc_out.data:
        fixed.append(inst.operation,
                     [qc_out.find_bit(q).index for q in inst.qubits])
    d = [None] * nq
    for i, p in enumerate(perm):
        d[p] = i
    nsw = 0
    for w in range(nq):
        while d[w] != w:
            j = d[w]
            fixed.swap(w, j)
            nsw += 1
            d[w], d[j] = d[j], d[w]
            if nsw > 4 * nq:
                raise RuntimeError("perm fix loop")
    return fixed, nsw


# ------------------------------------------------------------------- tools
def run_compactq(inp):
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from compactq import optimize_search
    circ = from_qiskit(inp)
    t0 = time.perf_counter()
    out = optimize_search(circ)
    dt = time.perf_counter() - t0
    return to_qiskit(out), dt


def run_qiskit(inp):
    from qiskit import transpile
    best, best_t = None, float("inf")
    for _ in range(3):
        t0 = time.perf_counter()
        out = transpile(inp, basis_gates=OUR_BASIS, optimization_level=3,
                        seed_transpiler=42)
        best_t = min(best_t, time.perf_counter() - t0)
        best = out
    return best, best_t


def run_pytket(inp):
    import pytket
    from pytket.extensions.qiskit import qiskit_to_tk, tk_to_qiskit
    from pytket.passes import (AutoRebase, CommuteThroughMultis,
                               FullPeepholeOptimise, RemoveRedundancies,
                               SequencePass)
    from pytket.predicates import CompilationUnit
    OT = pytket.circuit.OpType
    ourset = {OT.H, OT.X, OT.Y, OT.Z, OT.S, OT.Sdg, OT.T, OT.Tdg,
              OT.Rx, OT.Ry, OT.Rz, OT.U1, OT.SX, OT.CX, OT.CZ, OT.SWAP}
    tk = qiskit_to_tk(inp)
    cu = CompilationUnit(tk)
    t0 = time.perf_counter()
    seq = SequencePass([CommuteThroughMultis(), RemoveRedundancies(),
                        FullPeepholeOptimise(), AutoRebase(ourset),
                        CommuteThroughMultis(), RemoveRedundancies()])
    seq.apply(cu)
    dt = time.perf_counter() - t0
    return tk_to_qiskit(cu.circuit), dt


def run_cirq(inp):
    """Cirq receives the same circuit as OpenQASM text (its own importer).

    If the shared basis contains gates Cirq's importer does not know, the
    input is re-expressed once at level 0 in its importable subset - the
    optimizer input stays unitary-identical either way.
    """
    import cirq
    from qiskit import qasm2, transpile
    try:
        text = qasm2.dumps(inp)
        c = circuit_from_qasm_text(text)
    except Exception:
        cirq_basis = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry",
                      "rz", "u", "cx", "cz", "swap", "cp"]
        inp2 = transpile(inp, basis_gates=cirq_basis, optimization_level=0,
                         seed_transpiler=42)
        text = qasm2.dumps(inp2)
        c = circuit_from_qasm_text(text)
    t0 = time.perf_counter()
    opt = cirq.optimize_for_target_gateset(c, gateset=cirq.CZTargetGateset())
    dt = time.perf_counter() - t0
    # metrics with the same ASAP-depth definition the other tools report
    ops = list(opt.all_operations())
    total = len(ops)
    two_q = sum(1 for op in ops if len(op.qubits) == 2)
    frontier = {}
    for op in ops:
        t = max(frontier.get(q, 0) for q in op.qubits) + 1
        for q in op.qubits:
            frontier[q] = t
    depth = max(frontier.values(), default=0)
    m = SimpleMetrics(total, two_q, depth)
    return c, opt, m, dt


class SimpleMetrics:
    def __init__(self, total, two_q, depth):
        self.total, self.two_q, self.depth = total, two_q, depth

    def __iter__(self):
        return iter((self.total, self.two_q, self.depth))


def circuit_from_qasm_text(text: str):
    from cirq.contrib.qasm_import import circuit_from_qasm
    return circuit_from_qasm(text)


def cirq_unitary_vs(c_opt, ref, nq):
    order = sorted(c_opt.all_qubits())[::-1]
    # cirq's qubit_order is big-endian (first qubit = most significant);
    # reversing maps it onto qiskit Operator's little-endian wire order
    if len(order) != nq:
        return None
    return fid(c_opt.unitary(qubit_order=order), ref)


# -------------------------------------------------------------------- main
def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    cirq_ok = True
    try:
        import cirq  # noqa: F401
    except Exception as e:
        # the cirq COLUMN is optional: a broken or policy-blocked cirq
        # install must not kill the whole suite (skips only that column)
        cirq_ok = False
        log(f"cirq unavailable - its column is skipped ({str(e)[:60]})")
    try:
        import mqt.bench  # noqa: F401
    except Exception:
        log("mqt-bench not installed - pip install mqt-bench cirq ply")
        return 2

    log("tool versions: " + _versions())
    rows = []
    for name, size in CANDIDATES:
        tag = f"{name}_n{size}"
        qc, reason = gen_unitary_mqt(name, size)
        if qc is None:
            log(f"skip {tag}: {reason}")
            continue
        try:
            inp = normalize(qc)
        except Exception as e:
            log(f"skip {tag}: normalize: {type(e).__name__}: {str(e)[:60]}")
            continue
        if inp.num_qubits == 0 or len(inp.data) == 0:
            log(f"skip {tag}: empty after normalization")
            continue
        # flatten onto a single q register: identical gate list and unitary,
        # but QASM consumers that order qubits per-register instead of by
        # global wire index (cirq importer, pytket converter) get a
        # well-defined wire mapping
        from qiskit import QuantumCircuit
        flat = QuantumCircuit(inp.num_qubits)
        for inst in inp.data:
            flat.append(inst.operation,
                        [inp.find_bit(q).index for q in inst.qubits])
        inp = flat
        nq = inp.num_qubits
        ref = None
        if nq <= REF_QUBITS:
            from qiskit.quantum_info import Operator
            ref = Operator(inp).data

        row = {"name": tag, "nq": nq,
               "in": metrics(inp)}
        for key, fn in (("qz", run_compactq), ("qk", run_qiskit),
                        ("pt", run_pytket)):
            try:
                out, dt = fn(inp)
                m = metrics(out)
                f = fid(Operator(out).data, ref) if ref is not None else None
                if key == "qk" and ref is not None:
                    fixed, _ = reify_perm(out)
                    f = fid(Operator(fixed).data, ref)
                row[key] = (*m, round(dt * 1000), f)
            except Exception as e:
                row[key] = None
                row[key + "_err"] = f"{type(e).__name__}: {str(e)[:60]}"
        # cirq (separate: returns cirq objects); column is optional
        if cirq_ok:
            try:
                _, c_opt, m, dt = run_cirq(inp)
                f = cirq_unitary_vs(c_opt, ref, nq) if ref is not None else None
                row["cq"] = (*m, round(dt * 1000), f)
            except Exception as e:
                row["cq"] = None
                row["cq_err"] = f"{type(e).__name__}: {str(e)[:60]}"
        else:
            row["cq"] = None
            row["cq_err"] = "cirq unavailable"

        rows.append(row)
        log(_fmt_row(row))

    _summary(rows)

    # machine-readable artifact: every claim in results/ is backed by a
    # committed raw file (v0.2.4 scorecard contract)
    import json as _json
    from provenance import environment
    doc = {
        "generated": _utcnow(),
        "protocol": ("MQT Bench algorithm-level circuits, trailing "
                     "measurement stripped, level-0 normalized to the "
                     "shared basis, single flattened register; metrics "
                     "total/2q/depth + wall ms + referee fidelity "
                     "(qiskit Operator, circuits <= 8q); compactq's "
                     "time includes its whole-circuit proof"),
        "tool_versions": _versions(),
        "environment": environment(extra_deps=("pytket", "cirq")),
        "circuits": len(rows),
        "records": rows,
    }
    out = REPO / "results" / "mqtbench.json"
    out.write_text(_json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    md = ["# MQT Bench four-way comparison (generated)", "",
          f"Generated {doc['generated']} | {len(rows)} circuits | "
          f"referee: qiskit Operator (<= 8q)", "",
          "| circuit | nq | in t/2q/d | compactq | qiskit L3 | pytket | cirq |",
          "|---|---:|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['name']} | {r['nq']} | {r['in'][0]}/{r['in'][1]}/"
                  f"{r['in'][2]} | {_fmt_cell(r, 'qz')} | {_fmt_cell(r, 'qk')}"
                  f" | {_fmt_cell(r, 'pt')} | {_fmt_cell(r, 'cq')} |")
    (REPO / "results" / "mqtbench.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8")
    log(f"wrote {out} + mqtbench.md")
    return 0


def _utcnow():
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


def _versions():
    import qiskit
    import compactq
    parts = [f"compactq {compactq.__version__}", f"qiskit {qiskit.__version__}"]
    try:
        import pytket
        parts.append(f"pytket {pytket.__version__}")
    except Exception:
        pass
    try:
        import cirq
        parts.append(f"cirq {cirq.__version__}")
    except Exception:
        pass
    try:
        import mqt.bench
        parts.append("mqt-bench " + getattr(mqt.bench, "__version__", "?"))
    except Exception:
        pass
    return ", ".join(parts)


def _fmt_row(row):
    return (f"{row['name']}: in {row['in'][0]}/{row['in'][1]}/{row['in'][2]} | "
            f"compactq {_fmt_cell(row, 'qz')} | qk {_fmt_cell(row, 'qk')} | "
            f"pt {_fmt_cell(row, 'pt')} | cirq {_fmt_cell(row, 'cq')}")


def _fmt_cell(row, key):
    v = row.get(key)
    if v is None:
        err = row.get(key + "_err", "")
        return f"FAILED ({err})" if err else "FAILED"
    t, q, d, ms, f = v
    fs = f" fid={f:.6f}" if f is not None else ""
    return f"{t}/{q}/{d} ({ms}ms){fs}"


def _summary(rows):
    n = len(rows)
    if not n:
        log("no circuits ran")
        return
    log(f"\n=== MQT Bench: {n} circuits ===")
    for a, b, label in (("qz", "qk", "compactq vs qiskit L3"),
                        ("qz", "pt", "compactq vs pytket"),
                        ("qz", "cq", "compactq vs cirq"),
                        ("qk", "cq", "qiskit vs cirq")):
        wins = ties = losses = ineq = errs = 0
        for r in rows:
            va, vb = r.get(a), r.get(b)
            if va is None or vb is None:
                errs += 1
                continue
            if vb[4] is not None and vb[4] < 1 - 1e-6:
                ineq += 1
                continue  # don't count invalid competitor outputs as wins
            if va[1] < vb[1]:
                wins += 1
            elif va[1] == vb[1]:
                ties += 1
            else:
                losses += 1
        log(f"{label}: 2q win {wins}, tie {ties}, lose {losses}"
            + (f", competitor-inequivalent {ineq}" if ineq else "")
            + (f", tool-errors {errs}" if errs else ""))
    ineq_by_tool = {}
    for r in rows:
        for k in ("qz", "qk", "pt", "cq"):
            v = r.get(k)
            if v is not None and v[4] is not None and v[4] < 1 - 1e-6:
                ineq_by_tool.setdefault(k, []).append(r["name"])
    if ineq_by_tool:
        log("inequivalent outputs: " + "; ".join(
            f"{k}: {len(v)} ({', '.join(v)})" for k, v in ineq_by_tool.items()))


if __name__ == "__main__":
    sys.exit(main())
