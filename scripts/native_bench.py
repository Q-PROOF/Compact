"""Native-basis benchmark: Compact vs Qiskit L3 on machine-native 2q gates.

Same protocol as realbench (level-0 normalized identical inputs), but every
output is translated to a hardware-native basis: ecr/cz/iswap + rz/sx/x.
Circuits <= 8q are refereed with qiskit Operator for both tools.
Usage: python scripts/native_bench.py [ecr|cz|iswap] [--max-qubits N]
"""
from __future__ import annotations
import re, sys, warnings
from pathlib import Path

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np
from qiskit import transpile
from qiskit.quantum_info import Operator
from realbench import prepare_unitary, normalize_qiskit, normalize_harder, metrics


def fid(u, v):
    return abs(np.sum(np.conj(u) * v)) / u.shape[0]


def run(native: str, max_q: int):
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from compactq import optimize_search
    from compactq.native import rebase

    root = REPO / "third_party" / "QASMBench" / "small"
    files = []
    for d in sorted(root.iterdir()):
        if d.is_dir():
            for qf in sorted(d.glob("*.qasm")):
                if "_transpiled" not in qf.name and "_out" not in qf.name:
                    files.append(qf)
    basis1q = ["rz", "sx", "x"]
    rows = []
    for f in files:
        name = f.stem
        text = f.read_text()
        if re.search(r"^\s*if\s*\(", text, re.M):
            continue
        prepared, reason = prepare_unitary(text)
        if prepared is None:
            continue
        m = re.search(r"_n(\d+)", f.name)
        nq = int(m.group(1)) if m else prepared.num_qubits
        if nq > max_q:
            continue
        try:
            inp = normalize_qiskit(prepared)
        except Exception:
            inp = normalize_harder(prepared)
        ref = Operator(inp).data
        # Compact
        circ = from_qiskit(inp)
        opt = optimize_search(circ)
        rb = rebase(opt, native, oneq="rz-sx-x")
        rb_qc = to_qiskit(rb)
        f_c = fid(Operator(rb_qc).data, ref)
        m_c = metrics(rb_qc)
        # Qiskit L3 same target
        qk = transpile(inp, basis_gates=[native] + basis1q,
                       optimization_level=3, seed_transpiler=42)
        f_q = fid(Operator(qk).data, ref)
        m_q = metrics(qk)
        rows.append((name, m_c, f_c, m_q, f_q))
        flag_c = "" if f_c > 1 - 1e-6 else " INEQ!"
        flag_q = "" if f_q > 1 - 1e-6 else " INEQ!"
        print(f"{name:22s} Compact {m_c[0]}/{m_c[1]}/{m_c[2]}{flag_c} | "
              f"qk-L3 {m_q[0]}/{m_q[1]}/{m_q[2]}{flag_q}")
    # inequivalent competitor outputs are excluded from win/loss tallies
    valid = [r for r in rows if r[2] > 1 - 1e-6 and r[4] > 1 - 1e-6]
    ineq = len(rows) - len(valid)
    rows = valid
    if ineq:
        print(f"({ineq} outputs excluded from tallies: inequivalent under referee)")
    cw = sum(1 for r in rows if r[1][1] < r[3][1])
    ct = sum(1 for r in rows if r[1][1] == r[3][1])
    cl = [r[0] for r in rows if r[1][1] > r[3][1]]
    dw = sum(1 for r in rows if r[1][2] < r[3][2])
    print(f"\n=== {native.upper()} target: {len(rows)} circuits | 2q: Compact wins "
          f"{cw}, tie {ct}, loses {len(cl)} {cl} | depth wins {dw} ===")


if __name__ == "__main__":
    native = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "ecr"
    mq = 8
    if "--max-qubits" in sys.argv:
        mq = int(sys.argv[sys.argv.index("--max-qubits") + 1])
    run(native, mq)
