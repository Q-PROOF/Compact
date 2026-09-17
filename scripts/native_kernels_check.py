"""Timing comparison: native vs Python block unitary + end-to-end deep run."""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import compactq_native
from compactq import Circuit, Gate
from compactq.equivalence import unitary
from compactq.qiskit_bridge import from_qiskit
from compactq.optimize import optimize_deep
from qiskit import qasm2, transpile

rng = random = __import__("random").Random(5)


def flatten(ops):
    names, pf, pc, qf, qc = [], [], [], [], []
    for g in ops:
        names.append(g.name)
        pf.extend(g.params)
        pc.append(len(g.params))
        qf.extend(g.qubits)
        qc.append(len(g.qubits))
    return names, pf, pc, qf, qc


def rand_circ(n, k):
    ops = []
    for _ in range(k):
        r = rng.random()
        if r < 0.3:
            ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
        elif r < 0.4:
            ops.append(Gate("cp", (rng.uniform(0, 6),), tuple(rng.sample(range(n), 2))))
        else:
            name = rng.choice(["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx",
                               "sxdg", "rx", "ry", "rz", "p", "u3"])
            nparams = 3 if name == "u3" else 1 if name in ("rx", "ry", "rz", "p") else 0
            ops.append(Gate(name, tuple(rng.uniform(0, 6.28) for _ in range(nparams)),
                            (rng.randrange(n),)))
    return Circuit(n, ops)


circuits = [rand_circ(2, rng.randint(10, 40)) for _ in range(200)]

t0 = time.perf_counter()
for c in circuits:
    unitary(c)
t_py = time.perf_counter() - t0

t0 = time.perf_counter()
for c in circuits:
    args = flatten(c.ops)
    compactq_native.block_unitary(c.num_qubits, *args)
t_nat = time.perf_counter() - t0

print(f"200x random 2q blocks: python {t_py*1000:.0f} ms | native {t_nat*1000:.0f} ms "
      f"({t_py / t_nat:.1f}x)")

OUR_BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p",
             "sx", "u", "u3", "u1", "u2", "cx", "cz", "swap", "cp"]
raw = qasm2.loads((ROOT / "third_party/QASMBench/small/qaoa_n6/qaoa_n6.qasm").read_text(),
                  custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
inp = from_qiskit(transpile(raw.decompose(reps=3), basis_gates=OUR_BASIS,
                            optimization_level=0, seed_transpiler=42))
t0 = time.perf_counter()
res = optimize_deep(inp, verify=False)
dt = time.perf_counter() - t0
print(f"optimize_deep qaoa_n6 with native: {res.stats()} in {dt*1000:.0f} ms")
