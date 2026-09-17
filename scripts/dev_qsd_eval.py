"""3q QSD prototype v2: CSD via 4x4 SVD of the Q11 block, completed to the
full cosine-sine decomposition; validated against Qiskit synthesis costs."""
import sys, math, warnings
sys.path.insert(0, r"C:\Users\patha\.zcode\workspace\default\compactq")
sys.path.insert(0, r"C:\Users\patha\.zcode\workspace\default\compactq\scripts")
warnings.filterwarnings("ignore")
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import random_unitary
from compactq.qiskit_bridge import from_qiskit, to_qiskit
from compactq import optimize_search

def jacobi_svd_4x4(M, iters=100):
    """One-sided Jacobi SVD of a square complex matrix. Returns U_svd (4x4
    columns = right singular vectors), s (singular values), V_left (4x4)."""
    A = np.array(M, dtype=complex)
    n = A.shape[0]
    V = np.eye(n, dtype=complex)
    for _ in range(iters):
        off = 0.0
        for p in range(n - 1):
            for q_ in range(p + 1, n):
                app = np.vdot(A[:, p], A[:, p]).real
                aqq = np.vdot(A[:, q_], A[:, q_]).real
                apq = np.vdot(A[:, p], A[:, q_])
                if abs(apq) < 1e-15:
                    continue
                off += abs(apq)
                delta = (aqq - app) / (2.0 * apq)
                t = (1.0 if delta >= 0 else -1.0) / (abs(delta) + math.sqrt(1 + delta * delta))
                cth = 1.0 / math.sqrt(1 + t * t)
                sth = cth * t
                phase = apq / abs(apq)
                G = np.eye(n, dtype=complex)
                G[p, p] = G[q_, q_] = cth
                G[p, q_] = sth * np.conj(phase)
                G[q_, p] = -sth * phase
                A = A @ G
                V = V @ G
        if off < 1e-13:
            break
    s = np.array([np.linalg.norm(A[:, j]) for j in range(n)])
    order = np.argsort(-s)
    s = s[order]
    V = V[:, order]
    left = A[:, order] @ np.diag(1.0 / np.where(s > 1e-13, s, 1.0))
    return left, s, V

def csd_3q(U):
    """U: 8x8, split q0. Returns (L1,L2,c,s,R1,R2, phase_ok)."""
    Q11, Q12 = U[:4, :4], U[:4, 4:]
    Q21, Q22 = U[4:, :4], U[4:, 4:]
    L1, c, R1h = np.linalg.svd(Q11)[:3]   # Q11 = L1 diag(c) R1h
    R1 = R1h.conj().T
    s = np.sqrt(np.clip(1 - c ** 2, 0, 1))
    # Q21 = L2 diag(s) R1^dag  -> L2 = Q21 R1 diag(1/s)
    L2 = Q21 @ R1 @ np.diag(1.0 / np.where(s > 1e-13, s, 1.0))
    # Q12 = L1 diag(s) R2^dag -> R2 = Q12^dag L1 diag(1/s)  (then conj)
    R2 = (np.diag(1.0 / np.where(s > 1e-13, s, 1.0)) @ L1.conj().T @ Q12).conj().T
    # validate Q22 = L2 diag(c) R2^dag
    rec = L2 @ np.diag(c) @ R2.conj().T
    ok = np.abs(rec - Q22).max() < 1e-8
    return L1, L2, c, s, R1, R2, ok

rng = np.random.default_rng(11)
print("generic-3q validation:")
worst = 0.0
for trial in range(5):
    U = random_unitary(8, seed=int(rng.integers(1e6))).data
    L1, L2, c, s, R1, R2, ok = csd_3q(U)
    worst = max(worst, np.abs(c**2 + s**2 - 1).max())
    print(f"  trial {trial}: csd-complete={ok}, cos={np.round(c,3)}")
print(f"  worst c^2+s^2-1: {worst:.2e}")

# Qiskit reference: what does the best generic 3q synthesis cost?
U = random_unitary(8, seed=42)
qc = QuantumCircuit(3); qc.unitary(U, range(3))
t = transpile(qc, basis_gates=["u3", "cx"], optimization_level=3, seed_transpiler=3)
print(f"qiskit generic-3q cx count: {t.count_ops().get('cx', 0)} (optimal known ~14)")

# structured-window empirical check: does ANY suite circuit improve if we
# hand its full unitary to qiskit's best 3q/4q synthesizer? (only <=4q)
from realbench import prepare_unitary, normalize_qiskit
from pathlib import Path
f = Path(r"C:\Users\patha\.zcode\workspace\default\compactq\third_party\QASMBench\small\basis_trotter_n4\basis_trotter_n4.qasm")
prepared, _ = prepare_unitary(f.read_text())
inp = normalize_qiskit(prepared)
opt = optimize_search(from_qiskit(inp))
print(f"basis_trotter_n4: compactq 2q = {opt.two_qubit_count()} (12 cx-equiv with swaps)")
# qiskit L3 with u3,cx (its own best resynthesis of the SAME input):
t4 = transpile(inp, basis_gates=["u3", "cx"], optimization_level=3, seed_transpiler=5)
print(f"basis_trotter_n4: qiskit-L3(u3,cx) 2q = {t4.count_ops().get('cx', 0)} (+ elided perms)")
