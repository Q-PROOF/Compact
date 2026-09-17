"""Real-world gauntlet: realistic algorithm circuits through every compactq
entry point, with Qiskit Operator as an INDEPENDENT correctness referee.

Families (built in Qiskit, normalized to compactq's basis at level 0 - same
protocol as realbench):  QFT / inverse QFT, GHZ, W-state, Grover, QAOA
MaxCut rings, Heisenberg trotter blocks, phase-estimate cp ladders,
MCX/MCP controlled gates (raw, to exercise the expander), random
Cliffords, hardware-style ansaetze (RealAmplitudes / ZZFeatureMap),
quantum-volume-style random circuits, plus the QASMBench small suite fed
through compactq's own OpenQASM importer.

Entry points exercised per family:
  optimize / optimize_search / optimize_search(verify=False) /
  optimize_deep / approximate (contract check) / optimize_for (exact
  levels + default mixed levels) / approximate_for_target (budget=0
  exactness + free-budget sanity) / route_aware on line and heavy-hex
  couplings (edge-legality + exactness when restored) / expand_circuit /
  OpenQASM 2 + 3 round-trips / CLI subprocess / qiskit_bridge.compactq_pass /
  determinism (two runs -> identical op lists).

Referee: phase-insensitive average gate fidelity |Tr(A^dag B)|/d computed
from qiskit.quantum_info.Operator, never from compactq's own unitary code.

Exit code 0 iff every check passes.
"""
from __future__ import annotations

import io
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from qiskit import QuantumCircuit, transpile  # noqa: E402
from qiskit.quantum_info import Operator  # noqa: E402

from compactq import (approximate, optimize, optimize_deep, optimize_search,  # noqa: E402
                  from_qasm, from_qasm3, to_qasm, to_qasm3)
from compactq.mcx import expand_circuit  # noqa: E402
from compactq.qiskit_bridge import from_qiskit, compactq_pass, to_qiskit  # noqa: E402
from compactq.target import Target, approximate_for_target, optimize_for  # noqa: E402
from compactq.hardware import route_aware  # noqa: E402

OUR_BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p",
             "sx", "cx", "cz", "swap"]

EXACT_TOL = 1e-7


# ----------------------------------------------------------------- referee
def fid_of(u: np.ndarray, v: np.ndarray) -> float:
    """Phase-insensitive average gate fidelity |Tr(A^dag B)| / d."""
    return abs(np.sum(np.conj(u) * v)) / u.shape[0]


def op_of(qc) -> np.ndarray:
    return Operator(qc).data


# ------------------------------------------------------------------ inputs
def normalize(qc):
    out = transpile(qc, basis_gates=OUR_BASIS, optimization_level=0,
                    seed_transpiler=42)
    names = {inst.operation.name for inst in out.data}
    if not names <= set(OUR_BASIS) | {"barrier"}:
        out = transpile(qc.decompose(reps=3), basis_gates=OUR_BASIS,
                        optimization_level=0, seed_transpiler=42)
    return out


def build_families(max_n: int = 8):
    """(name, normalized_qiskit_circuit, referee_unitary) tuples."""
    fam = []
    for n in range(3, max_n + 1):
        qc = QuantumCircuit(n)
        for j in range(n):
            qc.h(j)
        for c in range(n - 1):
            for t in range(c + 1, n):
                qc.cp(np.pi / 2 ** (t - c), c, t)
        qc.swap(0, n - 1)
        fam.append((f"qft{n}", qc))
    for n in range(3, 7):
        qc = QuantumCircuit(n)
        for j in range(n):
            qc.h(j)
        for c in reversed(range(n - 1)):
            for t in reversed(range(c + 1, n)):
                qc.cp(-np.pi / 2 ** (t - c), c, t)
        qc.swap(0, n - 1)
        for j in range(n):
            qc.h(j)
        fam.append((f"iqft{n}", qc))
    for n in range(4, max_n + 1):
        qc = QuantumCircuit(n)
        qc.h(0)
        for j in range(1, n):
            qc.cx(j - 1, j)
        fam.append((f"ghz{n}", qc))
    for n in (4, 6):
        qc = QuantumCircuit(n)
        qc.ry(np.arccos(np.sqrt(1 / n)), 0)
        for k in range(1, n):
            qc.ry(np.arccos(np.sqrt(1 / (n - k + 1))), k)
            for j in range(k):
                qc.cz(k - 1, k)
        fam.append((f"wstate{n}", qc))
    for n in (3, 4, 5):
        qc = QuantumCircuit(n)
        qc.h(range(n))
        for it in range(2):
            qc.cz(0, n - 1)  # oracle-ish phase kickback
            for c in range(1, n - 1):
                qc.cx(c, n - 1)
            qc.h(range(n))
            for j in range(n):  # diffusion
                qc.rz(np.pi, j)
            qc.h(range(n))
            qc.x(range(n))
            qc.h(n - 1)
            qc.mcx(list(range(n - 1)), n - 1)
            qc.h(n - 1)
            qc.x(range(n))
            qc.h(range(n))
        fam.append((f"grover{n}", qc))
    for n in range(4, max_n + 1):
        for p in (1, 2, 3):
            qc = QuantumCircuit(n)
            rng = np.random.default_rng(100 + n)
            for j in range(n):
                qc.h(j)
            for _ in range(p):
                for j in range(n):
                    qc.rzz(2 * np.pi / 4, j, (j + 1) % n)
                    qc.rz(2 * rng.uniform(0.1, 0.5), j)
            fam.append((f"qaoa{n}_p{p}", qc))
    for n in (4, 6, 8):
        qc = QuantumCircuit(n)
        for layer in range(3):
            for j in range(n):
                qc.u(np.pi / 2, np.pi / 3, np.pi / 7, j)
            for j in range(0, n - 1, 2):
                qc.rzz(0.7, j, j + 1)
                qc.rzx(0.4, j, j + 1)
            for j in range(1, n - 1, 2):
                qc.rzz(0.7, j, j + 1)
                qc.rzx(0.4, j, j + 1)
        fam.append((f"trotter{n}", qc))
    for n in (4, 6):
        qc = QuantumCircuit(n)
        for j in range(n):
            qc.h(j)
        for sh in range(n - 1):
            qc.cp(np.pi / 2 ** (sh + 1), 0, n - 1 - sh)
        for j in range(n - 1):
            qc.cp(np.pi / 3, j, j + 1)
        fam.append((f"pwrap{n}", qc))
    for n in (4, 5, 6):
        qc = QuantumCircuit(n)
        qc.x(range(n - 1))
        qc.h(n - 1)
        qc.mcx(list(range(n - 1)), n - 1)
        qc.mcp(np.pi / 4, list(range(n - 1)), n - 1)
        qc.x(range(n - 1))
        fam.append((f"mcx{n}", qc))
    for n in (4, 6, 8):
        rng = np.random.default_rng(7 * n)
        from qiskit.quantum_info.random import random_clifford
        qc = random_clifford(n, seed=int(rng.integers(1, 10_000))).to_circuit()
        fam.append((f"cliff{n}", qc))
    for n in range(4, max_n + 1):
        for reps in (2, 3):
            qc = QuantumCircuit(n)
            rng = np.random.default_rng(1000 * n + reps)
            for r in range(reps):
                for j in range(n):
                    qc.ry(rng.uniform(0, np.pi), j)
                    qc.rz(rng.uniform(0, 2 * np.pi), j)
                for j in range(n - 1):
                    qc.cx(j, j + 1)
            for j in range(n):
                qc.ry(rng.uniform(0, np.pi), j)
            fam.append((f"ansatz{n}_r{reps}", qc))
    for n in (4, 6, 8):
        qc = QuantumCircuit(n)
        rng = np.random.default_rng(31)
        for j in range(n):
            qc.h(j)
        for j in range(n - 1):
            qc.cp(2 * rng.uniform(0.1, 1.0), j, j + 1)
        for r in range(2):
            for j in range(n):
                qc.u(rng.uniform(0, np.pi), 0, 0, j)
            for j in range(n - 1):
                qc.cp(np.pi / 2 ** (r + 1), j, j + 1)
        fam.append((f"zzfeat{n}", qc))
    for n in (6, 8):
        qc = QuantumCircuit(n)
        rng = np.random.default_rng(5)
        for layer in range(n):
            perm = rng.permutation(n)
            for j in range(0, n - 1, 2):
                a, b = int(perm[j]), int(perm[j + 1])
                qc.u(rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi), 0, a)
                qc.u(rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi), 0, b)
                qc.cz(a, b)
        fam.append((f"qvol{n}", qc))
    out = []
    for name, qc in fam:
        norm = normalize(qc)
        out.append((name, norm, op_of(norm)))
    return out


# -------------------------------------------------------------- validation
def validate(c: "object", tag: str, errs: list):
    from compactq import Circuit
    if not isinstance(c, Circuit):
        errs.append(f"{tag}: not a compactq Circuit ({type(c)})")
        return
    for i, g in enumerate(c.ops):
        for q in g.qubits:
            if not (0 <= q < c.num_qubits):
                errs.append(f"{tag}: op {i} {g.name} qubit {q} out of range")
        for p in g.params:
            if not np.isfinite(p):
                errs.append(f"{tag}: op {i} {g.name} non-finite param {p}")
        if g.name not in KNOWN_GATES:
            errs.append(f"{tag}: op {i} unknown gate name {g.name!r}")


KNOWN_GATES = {"h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz",
               "p", "sx", "sxdg", "u3", "u", "cx", "cz", "swap", "cp", "cu1",
               "id", "rzz", "rzx", "mcx", "mcp", "barrier", "reset", "measure", "ecr", "iswap",
               "ch", "crz", "cu3", "csx", "ccx", "cy", "crx", "cry"}


def exact_check(name, entry, out_c, ref_u, res):
    errs = []
    validate(out_c, f"{name}/{entry}", errs)
    qc_out = to_qiskit(out_c)
    f = fid_of(op_of(qc_out), ref_u)
    if f < 1 - EXACT_TOL:
        errs.append(f"{name}/{entry}: EXACTNESS FAIL fid={f:.12f} "
                    f"({len(out_c.ops)} ops, {out_c.two_qubit_count()} 2q)")
    res.append((name, entry, "FAIL" if errs else "PASS", errs or
                (None, f"fid={f:.10f}")[1]))


# ------------------------------------------------------------------ runner
FAILS = []
ROWS = []


def run_case(name, qc, ref_u):
    circ = from_qiskit(qc)
    t2 = (circ.two_qubit_count(), len(circ.ops), circ.depth())

    def rec(entry, status, detail=""):
        ROWS.append((name, entry, status, detail))
        if status == "FAIL":
            for e in (detail if isinstance(detail, list) else [detail]):
                FAILS.append((name, entry, e))
            print(f"  FAIL {name}/{entry}: {detail}")

    # 1. optimize
    try:
        t0 = time.perf_counter()
        out = optimize(circ)
        dt = time.perf_counter() - t0
        errs = []
        validate(out, f"{name}/optimize", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        rec("optimize", "FAIL" if errs else "PASS",
            errs or f"{t2[0]}->{out.two_qubit_count()} 2q, fid={f:.9f}, {dt*1000:.0f}ms")
    except Exception as e:
        rec("optimize", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 2. optimize_search default
    try:
        t0 = time.perf_counter()
        out = optimize_search(circ)
        dt = time.perf_counter() - t0
        errs = []
        validate(out, f"{name}/search", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        rec("search", "FAIL" if errs else "PASS",
            errs or f"{t2[0]}->{out.two_qubit_count()} 2q, fid={f:.9f}, {dt*1000:.0f}ms")
    except Exception as e:
        rec("search", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 3. optimize_search(verify=False)  <- historically the dangerous path
    try:
        out = optimize_search(circ, verify=False)
        errs = []
        validate(out, f"{name}/search-nv", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        rec("search-nv", "FAIL" if errs else "PASS",
            errs or f"{t2[0]}->{out.two_qubit_count()} 2q, fid={f:.9f}")
    except Exception as e:
        rec("search-nv", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 4. optimize_deep
    try:
        out = optimize_deep(circ)
        errs = []
        validate(out, f"{name}/deep", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        rec("deep", "FAIL" if errs else "PASS",
            errs or f"{t2[0]}->{out.two_qubit_count()} 2q, fid={f:.9f}")
    except Exception as e:
        rec("deep", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 5. approximate: documented bound infid <= #2q_blocks * (1 - min_fid)
    try:
        mf = 0.99
        out = approximate(circ, min_fidelity=mf)
        errs = []
        validate(out, f"{name}/approx", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        bound = 1 - out.two_qubit_count() * (1 - mf)
        if f < bound - 1e-6:
            errs.append(f"fid={f:.8f} < documented bound {bound:.8f}")
        if f < 0.9:
            errs.append(f"fid={f:.6f} catastrophic")
        rec("approx", "FAIL" if errs else "PASS",
            errs or f"{t2[0]}->{out.two_qubit_count()} 2q, fid={f:.6f} >= bound {bound:.4f}")
    except Exception as e:
        rec("approx", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 6a. optimize_for EXACT levels only -> must be exactly equivalent
    try:
        tgt = Target(cx_fidelity={(0, 1): 0.995, (1, 0): 0.98,
                                  frozenset((1, 2)): 0.993},
                     default_cx_fidelity=0.99)
        out, est = optimize_for(circ, tgt, approx_levels=(1.0,))
        errs = []
        validate(out, f"{name}/for-exact", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        rec("for-exact", "FAIL" if errs else "PASS",
            errs or f"est={est:.6f} fid={f:.9f}")
    except Exception as e:
        rec("for-exact", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 6b. optimize_for default mixed levels: model must not overpromise
    try:
        tgt = Target(default_cx_fidelity=0.99)
        out, est = optimize_for(circ, tgt)
        errs = []
        validate(out, f"{name}/for-mixed", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < est - 0.02:
            errs.append(f"measured fid={f:.6f} << model est={est:.6f}")
        rec("for-mixed", "FAIL" if errs else "PASS",
            errs or f"est={est:.6f} measured={f:.6f}")
    except Exception as e:
        rec("for-mixed", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 7a. approximate_for_target budget=0 -> exact
    try:
        tgt = Target(default_cx_fidelity=0.99)
        out, est = approximate_for_target(circ, tgt, budget=0.0)
        errs = []
        validate(out, f"{name}/aft0", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"budget=0 not exact: fid={f:.12f}")
        rec("aft0", "FAIL" if errs else "PASS", errs or f"fid={f:.9f}")
    except Exception as e:
        rec("aft0", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 7b. approximate_for_target free budget: measured >= est - 0.02
    try:
        tgt = Target(default_cx_fidelity=0.99)
        out, est = approximate_for_target(circ, tgt)
        errs = []
        validate(out, f"{name}/aft", errs)
        f = fid_of(op_of(to_qiskit(out)), ref_u)
        if f < est - 0.02:
            errs.append(f"measured fid={f:.6f} << est={est:.6f}")
        rec("aft", "FAIL" if errs else "PASS",
            errs or f"est={est:.6f} measured={f:.6f}")
    except Exception as e:
        rec("aft", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 8. route_aware: line + heavy-hex-ish coupling
    n = circ.num_qubits
    line = [(j, j + 1) for j in range(n - 1)]
    hexish = line + [(j, j + 2) for j in range(0, n - 2, 4)]
    for tag, coup in (("route-line", line), ("route-hex", hexish)):
        try:
            routed, pos = route_aware(circ, coup, restore=True)
            errs = []
            validate(routed, f"{name}/{tag}", errs)
            edges = {frozenset(e) for e in coup}
            for g in routed.ops:
                if len(g.qubits) == 2 and frozenset(g.qubits) not in edges:
                    errs.append(f"{g.name}{g.qubits} off-coupling")
                    break
            if sorted(pos.keys()) != list(range(n)) or sorted(pos.values()) != list(range(n)):
                errs.append(f"final_position not a permutation: {pos}")
            f = fid_of(op_of(to_qiskit(routed)), ref_u)
            if f < 1 - EXACT_TOL:
                errs.append(f"fid={f:.12f}")
            rec(tag, "FAIL" if errs else "PASS",
                errs or f"{t2[0]}->{routed.two_qubit_count()} 2q, fid={f:.9f}")
        except Exception as e:
            rec(tag, "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 9. qasm2 round trip (text -> compactq -> text -> qiskit)
    try:
        errs = []
        text = to_qasm(circ)
        back = from_qasm(text)
        f = fid_of(op_of(to_qiskit(back)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"qasm2 round-trip fid={f:.12f}")
        rec("qasm2-rt", "FAIL" if errs else "PASS", errs or "ok")
    except Exception as e:
        rec("qasm2-rt", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 10. qasm3 round trip
    try:
        errs = []
        text = to_qasm3(circ)
        back = from_qasm3(text)
        f = fid_of(op_of(to_qiskit(back)), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"qasm3 round-trip fid={f:.12f}")
        rec("qasm3-rt", "FAIL" if errs else "PASS", errs or "ok")
    except Exception as e:
        rec("qasm3-rt", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 11. CLI subprocess
    try:
        errs = []
        with tempfile.TemporaryDirectory() as td:
            inp = Path(td) / "in.qasm"
            outp = Path(td) / "out.qasm"
            inp.write_text(to_qasm(circ))
            r = subprocess.run(
                [sys.executable, "-m", "compactq", str(inp), "-o", str(outp)],
                capture_output=True, text=True, cwd=str(REPO), timeout=600)
            if r.returncode != 0:
                errs.append(f"CLI rc={r.returncode}: {r.stderr[-400:]}")
            elif not outp.exists():
                errs.append("CLI produced no output file")
            else:
                back = from_qasm(outp.read_text())
                f = fid_of(op_of(to_qiskit(back)), ref_u)
                if f < 1 - EXACT_TOL:
                    errs.append(f"CLI output fid={f:.12f}")
        rec("cli", "FAIL" if errs else "PASS", errs or "ok")
    except Exception as e:
        rec("cli", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 12. qiskit plugin pass
    try:
        errs = []
        out_qc = compactq_pass(qc)
        f = fid_of(op_of(out_qc), ref_u)
        if f < 1 - EXACT_TOL:
            errs.append(f"plugin fid={f:.12f}")
        rec("plugin", "FAIL" if errs else "PASS", errs or "ok")
    except Exception as e:
        rec("plugin", "FAIL", f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    # 13. determinism
    try:
        a = optimize_search(circ)
        b = optimize_search(circ)
        same = [(g.name, g.params, g.qubits) for g in a.ops] == \
               [(g.name, g.params, g.qubits) for g in b.ops]
        rec("determinism", "PASS" if same else "FAIL",
            "" if same else "two runs differ")
    except Exception as e:
        rec("determinism", "FAIL", f"{type(e).__name__}: {e}")


def run_mcx_raw():
    """MCX family through the compactq expander (no qiskit normalization)."""
    for n in (4, 5, 6):
        qc = QuantumCircuit(n)
        qc.x(range(n - 1))
        qc.h(n - 1)
        qc.mcx(list(range(n - 1)), n - 1)
        qc.mcp(np.pi / 4, list(range(n - 1)), n - 1)
        qc.x(range(n - 1))
        ref = op_of(qc)
        circ = from_qiskit(qc)
        expanded = expand_circuit(circ)
        out = optimize_search(expanded)
        f = fid_of(op_of(to_qiskit(out)), ref)
        ok = f >= 1 - EXACT_TOL
        ROWS.append((f"mcxraw{n}", "expand+search", "PASS" if ok else "FAIL",
                     f"fid={f:.9f}"))
        if not ok:
            FAILS.append((f"mcxraw{n}", "expand+search", f"fid={f:.12f}"))
        print(f"  mcxraw{n}: fid={f:.9f} "
              f"({circ.two_qubit_count()} raw 2q -> {out.two_qubit_count()} 2q)")


def run_qasmbench_importer():
    """QASMBench small suite through compactq's OWN OpenQASM importer."""
    root = REPO / "third_party" / "QASMBench" / "small"
    from realbench import prepare_unitary  # reuse skip logic
    import re
    files = []
    for d in sorted(root.iterdir()):
        if d.is_dir():
            for qf in d.glob("*.qasm"):
                if "_transpiled" not in qf.name and "_out" not in qf.name:
                    files.append(qf)
                    break
    for f in files:
        name = f.parent.name
        text = f.read_text()
        if re.search(r"^\s*if\s*\(", text, re.M):
            continue
        prepared, reason = prepare_unitary(text)
        if prepared is None:
            continue
        m = re.search(r"_n(\d+)", f.name)
        nq = int(m.group(1)) if m else prepared.num_qubits
        if nq > 8:
            continue
        ref = op_of(prepared)
        try:
            circ = from_qasm(text)
            if circ.num_qubits != prepared.num_qubits:
                # importer may size the register differently; pad referee
                pass
            out = optimize_search(circ)
            errs = []
            validate(out, f"qb:{name}/import+search", errs)
            f = fid_of(op_of(to_qiskit(out)), ref)
            if f < 1 - EXACT_TOL:
                errs.append(f"fid={f:.12f}")
            ok = not errs
            ROWS.append((f"qb:{name}", "import+search",
                         "PASS" if ok else "FAIL", errs or
                         f"{circ.two_qubit_count()}->{out.two_qubit_count()} 2q fid={f:.9f}"))
            if not ok:
                FAILS.append((f"qb:{name}", "import+search", str(errs)))
            print(f"  qb:{name}: fid={f:.9f} "
                  f"({circ.two_qubit_count()}->{out.two_qubit_count()} 2q)")
        except Exception as e:
            FAILS.append((f"qb:{name}", "import+search",
                          f"{type(e).__name__}: {e}\n{traceback.format_exc()}"))
            ROWS.append((f"qb:{name}", "import+search", "FAIL",
                         f"{type(e).__name__}: {e}"))
            print(f"  qb:{name}: EXCEPTION {type(e).__name__}: {e}")


def run_big_no_verify():
    """>8q circuits: verify-off path must still be SOUND (10q referee)."""
    for n in (9, 10):
        qc = QuantumCircuit(n)
        for j in range(n):
            qc.h(j)
        for c in range(n - 1):
            for t in range(c + 1, n):
                qc.cp(np.pi / 2 ** (t - c), c, t)
        for layer in range(2):
            for j in range(n - 1):
                qc.rzz(0.6, j, j + 1)
                qc.ry(0.3, j)
        norm = normalize(qc)
        ref = op_of(norm)
        circ = from_qiskit(norm)
        out = optimize_search(circ, verify=False)
        errs = []
        validate(out, f"big{n}/search-nv", errs)
        f = fid_of(op_of(to_qiskit(out)), ref)
        if f < 1 - EXACT_TOL:
            errs.append(f"fid={f:.12f}")
        ok = not errs
        ROWS.append((f"big{n}", "search-nv", "PASS" if ok else "FAIL",
                     errs or f"{circ.two_qubit_count()}->{out.two_qubit_count()} 2q fid={f:.9f}"))
        if not ok:
            FAILS.append((f"big{n}", "search-nv", str(errs)))
        print(f"  big{n}: fid={f:.9f} "
              f"({circ.two_qubit_count()}->{out.two_qubit_count()} 2q)")
        # CLI dispatch on the same width: the entry point must route wide
        # circuits to randomized verification instead of crashing (external
        # audit found 18 CLI failures here that in-process checks missed)
        try:
            import tempfile
            r, fcli, cli_ok = None, -1.0, False
            with tempfile.TemporaryDirectory() as td:
                inp = Path(td) / "in.qasm"
                outp = Path(td) / "out.qasm"
                inp.write_text(to_qasm(circ))
                r = subprocess.run(
                    [sys.executable, "-m", "compactq", str(inp), "-o", str(outp)],
                    capture_output=True, text=True, cwd=str(REPO), timeout=900)
                cli_ok = (r.returncode == 0 and outp.exists()
                          and "proof status" in r.stderr)
            if cli_ok and outp.exists():
                back = from_qasm(outp.read_text())
                fcli = fid_of(op_of(to_qiskit(back)), ref)
                cli_ok = fcli > 1 - 1e-6
            if cli_ok:
                detail = f"rc={r.returncode}, fid={fcli:.9f}"
            else:
                detail = f"rc={r.returncode}: {r.stderr[-200:] if r else 'no run'}"
            ROWS.append((f"big{n}", "cli-dispatch",
                         "PASS" if cli_ok else "FAIL", detail))
            if not cli_ok:
                FAILS.append((f"big{n}", "cli-dispatch", r.stderr[-300:]))
            verdict = "PASS" if cli_ok else "FAIL"
            print(f"  big{n}/cli: {verdict}")
        except Exception as e:
            FAILS.append((f"big{n}", "cli-dispatch", f"{type(e).__name__}: {e}"))
            ROWS.append((f"big{n}", "cli-dispatch", "FAIL", str(e)))
            print(f"  big{n}/cli: EXCEPTION {type(e).__name__}")


def run_edge_cases():
    res = []

    def check(label, fn):
        try:
            fn()
            res.append((label, True, ""))
        except Exception as e:
            res.append((label, False, f"{type(e).__name__}: {e}\n{traceback.format_exc()}"))

    def empty():
        qc = QuantumCircuit(2)
        out = optimize_search(from_qiskit(qc))
        assert len(out.ops) == 0

    def one_qubit():
        qc = QuantumCircuit(1)
        qc.rz(0.3, 0)
        qc.ry(1.1, 0)
        qc.rz(0.9, 0)
        qc.h(0)
        qc.t(0)
        ref = op_of(qc)
        out = optimize_search(from_qiskit(qc))
        f = fid_of(op_of(to_qiskit(out)), ref)
        assert f >= 1 - EXACT_TOL, f"fid={f}"

    def weird_names():
        text = """
OPENQASM 2.0;
include "qelib1.inc";
qreg data[3];
qreg anc[1];
creg c[3];
h data[0];
cx data[0], data[1];
cx data[1], data[2];
cz anc[0], data[2];
t data;
measure data -> c;
"""
        circ = from_qasm(text)
        assert circ.num_qubits == 4, circ.num_qubits
        out = optimize(circ)
        assert out.num_qubits == 4

    def qasm3_import():
        text = """
OPENQASM 3.0;
include "stdgates.inc";
qubit[4] q;
bit[4] c;
h q;
cx q[0], q[1];
cp(pi/4) q[1], q[2];
cz q[2], q[3];
"""
        ref_qc = QuantumCircuit(4)
        ref_qc.h(range(4))
        ref_qc.cx(0, 1)
        ref_qc.cp(np.pi / 4, 1, 2)
        ref_qc.cz(2, 3)
        ref = op_of(ref_qc)
        circ = from_qasm3(text)
        out = optimize(circ)
        f = fid_of(op_of(to_qiskit(out)), ref)
        assert f >= 1 - EXACT_TOL, f"fid={f}"

    def degenerate_angles():
        qc = QuantumCircuit(2)
        qc.rz(0.0, 0)
        qc.rz(2 * np.pi, 1)
        qc.cx(0, 1)
        qc.cx(0, 1)
        qc.rx(0.0, 0)
        qc.ry(4 * np.pi, 1)
        ref = op_of(qc)
        out = optimize_search(from_qiskit(qc))
        f = fid_of(op_of(to_qiskit(out)), ref)
        assert f >= 1 - EXACT_TOL, f"fid={f}"

    def duplicate_ops():
        qc = QuantumCircuit(3)
        for _ in range(6):
            qc.cx(0, 1)
            qc.h(2)
        ref = op_of(qc)
        out = optimize_search(from_qiskit(qc))
        f = fid_of(op_of(to_qiskit(out)), ref)
        assert f >= 1 - EXACT_TOL, f"fid={f}"

    check("edge/empty", empty)
    check("edge/one-qubit", one_qubit)
    check("edge/weird-names", weird_names)
    check("edge/qasm3-import", qasm3_import)
    check("edge/degenerate-angles", degenerate_angles)
    check("edge/duplicate-ops", duplicate_ops)
    for label, ok, detail in res:
        ROWS.append((label, "-", "PASS" if ok else "FAIL", detail or "ok"))
        if not ok:
            FAILS.append((label, "-", detail))
            print(f"  FAIL {label}: {detail}")
        else:
            print(f"  {label}: ok")


# -------------------------------------------------------------------- main
def main():
    t_start = time.perf_counter()
    print("building families (qiskit normalization + referee unitaries)...")
    fam = build_families(max_n=8)
    print(f"{len(fam)} families; running gauntlet...\n")

    heavy = {"qft8", "qaoa8_p3", "trotter8", "ansatz8_r3", "cliff8",
             "zzfeat8", "qvol8", "qaoa7_p3", "ghz8", "ansatz7_r3",
             "qaoa7_p2", "qaoa8_p2", "ansatz7_r2", "qaoa6_p3"}
    for name, qc, ref in fam:
        t0 = time.perf_counter()
        print(f"[{name}] in={qc.size()}g/{qc.num_qubits}q")
        run_case(name, qc, ref)
        print(f"  ({time.perf_counter() - t0:.1f}s)")

    print("\n[MCX raw expander]")
    run_mcx_raw()
    print("\n[QASMBench through compactq's own importer]")
    run_qasmbench_importer()
    print("\n[>8q no-verify soundness]")
    run_big_no_verify()
    print("\n[edge cases]")
    run_edge_cases()

    dt = time.perf_counter() - t_start
    n_pass = sum(1 for r in ROWS if r[2] == "PASS")
    n_fail = sum(1 for r in ROWS if r[2] == "FAIL")
    print(f"\n=== GAUNTLET: {n_pass} PASS, {n_fail} FAIL, {len(ROWS)} checks, "
          f"{dt:.0f}s ===")
    if FAILS:
        print("\nFAILURES:")
        for name, entry, detail in FAILS:
            d = detail.splitlines()[0] if detail else ""
            print(f"  {name}/{entry}: {d[:300]}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
