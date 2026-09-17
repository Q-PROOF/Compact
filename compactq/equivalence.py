"""Exact equivalence checking (pure Python, for small qubit counts).

Builds the full unitary of a circuit and computes the Hilbert-Schmidt
fidelity |Tr(U1^dag U2)| / d^2, which equals 1 iff the two circuits are
equivalent up to global phase.
"""
from __future__ import annotations

import cmath
import math
import operator

from .circuit import Circuit, Gate
from .linalg import gate_matrix

try:  # optional Rust acceleration raises the proof ceiling to 8 qubits
    import compactq_native as _NATIVE
except Exception:
    _NATIVE = None

_MAX_QUBITS = 8 if _NATIVE is not None else 6


def _apply_1q(U, mat, q, dim):
    bit = 1 << q
    a, b, c, d = mat
    for r0 in range(dim):
        if r0 & bit:
            continue
        r1 = r0 | bit
        x0, x1 = U[r0], U[r1]
        U[r0] = [a * x0[k] + b * x1[k] for k in range(dim)]
        U[r1] = [c * x0[k] + d * x1[k] for k in range(dim)]


def _apply_cx(U, ctrl, tgt, dim):
    cb, tb = 1 << ctrl, 1 << tgt
    for r in range(dim):
        if (r & cb) and not (r & tb):
            U[r], U[r | tb] = U[r | tb], U[r]


def _apply_cz(U, q0, q1, dim):
    b0, b1 = 1 << q0, 1 << q1
    for r in range(dim):
        if (r & b0) and (r & b1):
            U[r] = [-x for x in U[r]]


def _apply_swap(U, q0, q1, dim):
    b0, b1 = 1 << q0, 1 << q1
    for r in range(dim):
        if (r & b0) and not (r & b1):
            U[r], U[r - b0 + b1] = U[r - b0 + b1], U[r]


def _flat_ops(circ: Circuit):
    names, pf, pc, qf, qc = [], [], [], [], []
    for g in circ.ops:
        names.append(g.name)
        pf.extend(g.params)
        pc.append(len(g.params))
        qf.extend(g.qubits)
        qc.append(len(g.qubits))
    return names, pf, pc, qf, qc


def _unitary_native(circ: Circuit):
    """Full unitary via the Rust kernels (sparse gate application)."""
    if _NATIVE is None:
        raise RuntimeError("native kernels unavailable")
    n = circ.num_qubits
    names, pf, pc, qf, qc = _flat_ops(circ)
    flat = _NATIVE.sim_unitary(n, names, pf, pc, qf, qc)
    dim = 1 << n
    return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1])
             for j in range(dim)]
            for i in range(dim)]


def unitary(circ: Circuit):
    """Full unitary of the circuit as a dim x dim list of row-lists.

    Uses the Rust kernel whenever it is available and supports every gate in
    the circuit (validated against the Python reference to ~1e-16); any
    unsupported gate falls back to the pure-Python reference implementation.
    """
    n = circ.num_qubits
    if n > _MAX_QUBITS:
        raise ValueError(f"unitary check limited to {_MAX_QUBITS} qubits, got {n}")
    if _NATIVE is not None:
        try:
            return _unitary_native(circ)
        except Exception:
            pass  # gate outside the kernel's table: pure-Python reference path
    dim = 1 << n
    U = [[1.0 if i == j else 0.0 for j in range(dim)] for i in range(dim)]
    for g in circ.ops:
        if len(g.qubits) == 1:
            _apply_1q(U, gate_matrix(g), g.qubits[0], dim)
        elif g.name == "cx":
            _apply_cx(U, g.qubits[0], g.qubits[1], dim)
        elif g.name == "cz":
            _apply_cz(U, g.qubits[0], g.qubits[1], dim)
        elif g.name == "swap":
            _apply_swap(U, g.qubits[0], g.qubits[1], dim)
        elif g.name in ("ecr", "iswap"):
            for op in _native_ops(g.name):
                if op.name == "cx":
                    _apply_cx(U, op.qubits[0], op.qubits[1], dim)
                else:
                    _apply_1q(U, gate_matrix(op), op.qubits[0], dim)
        elif g.name == "cp":
            c, t = g.qubits
            b0, b1 = 1 << c, 1 << t
            ph = cmath.exp(1j * g.params[0])
            for r in range(dim):
                if (r & b0) and (r & b1):
                    U[r] = [x * ph for x in U[r]]
        else:
            raise ValueError(f"no unitary rule for {g.name!r}")
    return U


def _sig(circ: Circuit):
    """Value signature of a circuit's op list (for the unitary memo)."""
    return tuple((g.name, g.params, g.qubits) for g in circ.ops)


_UNITARY_CACHE: dict = {}


def _unitary_native_flat(circ: Circuit):
    """Raw interleaved re/im unitary from the Rust kernel (no conversion)."""
    n = circ.num_qubits
    names, pf, pc, qf, qc = _flat_ops(circ)
    return _NATIVE.sim_unitary(n, names, pf, pc, qf, qc), n


def _fidelity(a_flat, b_flat, dim: int) -> float:
    """|Tr(A^dag B)| / dim over interleaved re/im flat unitaries.

    Uses the Rust trace2 kernel when available; the Python fallback keeps
    the four sums at C level via map/sum over strided slices.
    """
    if _NATIVE is not None and hasattr(_NATIVE, "trace2"):
        return _NATIVE.trace2(a_flat, b_flat)
    are = a_flat[0::2]
    aim = a_flat[1::2]
    bre = b_flat[0::2]
    bim = b_flat[1::2]
    mul = operator.mul
    re = sum(map(mul, are, bre)) + sum(map(mul, aim, bim))
    im = sum(map(mul, aim, bre)) - sum(map(mul, are, bim))
    return math.hypot(re, im) / dim


def _cached_flat_unitary(circ: Circuit):
    """Unitary (flat, native) with a small value-keyed memo; the whole-circuit
    reference unitary is recomputed for every candidate verification, so this
    is the hot path in the search."""
    if _NATIVE is None:
        return None
    key = (circ.num_qubits, _sig(circ))
    hit = _UNITARY_CACHE.get(key)
    if hit is not None:
        return hit
    try:
        flat, _ = _unitary_native_flat(circ)
    except Exception:
        return None  # unsupported gate: caller uses the reference path
    if len(_UNITARY_CACHE) > 16:
        _UNITARY_CACHE.clear()
    _UNITARY_CACHE[key] = flat
    return flat


def fidelity(circ_a: Circuit, circ_b: Circuit) -> float:
    """|Tr(A^dag B)| / d in [0, 1]; 1.0 == equivalent up to global phase."""
    if circ_a.num_qubits != circ_b.num_qubits:
        raise ValueError("qubit count mismatch")
    fa = _cached_flat_unitary(circ_a)
    if fa is not None:
        fb = _cached_flat_unitary(circ_b)
        if fb is not None:
            return _fidelity(fa, fb, 1 << circ_a.num_qubits)
    A = unitary(circ_a)
    B = unitary(circ_b)
    dim = len(A)
    tr = 0j
    for i in range(dim):
        row_a = A[i]
        row_b = B[i]
        for j in range(dim):
            tr += row_a[j].conjugate() * row_b[j]
    return abs(tr) / dim


def check_equivalent(circ_a: Circuit, circ_b: Circuit, tol: float = 1e-7) -> bool:
    return fidelity(circ_a, circ_b) > 1.0 - tol

def _native_ops(name):
    """{u3 + cx} expansion of a native gate on wires (0, 1)."""
    from .native import NATIVE_EXPANSION
    return NATIVE_EXPANSION[name]
