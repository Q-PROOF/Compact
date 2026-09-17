"""Quickstart: optimize a QFT circuit with compactq (optionally compare to Qiskit)."""
import compactq
from compactq.benchmarks import qft

circ = qft(4)
opt = compactq.optimize(circ)

print(f"raw : {circ.stats()}")
print(f"compactq: {opt.stats()}")

try:
    from qiskit import transpile, qasm2 as qk_qasm
    qc = qk_qasm.loads(compactq.to_qasm(circ))
    basis = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p", "cx", "cz", "swap"]
    out = transpile(qc, basis_gates=basis, optimization_level=3, seed_transpiler=42)
    cx = out.count_ops().get("cx", 0)
    print(f"qiskit L3: {int(sum(out.count_ops().values()))} gates, {cx} cx, depth {out.depth()}")
except ImportError:
    print("(install qiskit for the comparison line: pip install qiskit)")
