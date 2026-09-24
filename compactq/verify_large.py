"""Large-circuit exactness verification by randomized state application.

For unitaries A, B the operator difference D = A - B obeys: if
``||D|psi>|| = 0`` on K random product states spanning enough randomness,
then ||D|| is negligible with overwhelming probability (the standard
randomized-verification argument used in cycle benchmarking).  We apply
both circuits to K=32 random product states with an exact numpy
statevector simulator and require agreement to 1e-8 on every state.

This extends the proof net beyond the 8-qubit dense ceiling to ~30 qubits
on today's laptops, at which point memory becomes the binding constraint.

OPTIONAL accelerator: needs numpy.  compactq's core remains zero-dependency.
"""
from __future__ import annotations

import time

RAND_MAX_QUBITS = 30        # 2^30 * 16 B per statevector: memory ceiling


def _np():
    import numpy
    return numpy


def _mat2q(name, params):
    np = _np()
    if name == "cx":
        return np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
                        dtype=complex)
    if name == "cz":
        return np.diag([1, 1, 1, -1]).astype(complex)
    if name == "swap":
        return np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]],
                        dtype=complex)
    if name == "cp":
        return np.diag([1, 1, 1, np.exp(1j * float(params[0]))]).astype(complex)
    if name in ("ecr", "iswap"):
        from .native import NATIVE_EXPANSION
        from .linalg import gate_matrix
        M = np.eye(4, dtype=complex)

        def apply_col(v):
            for gg in NATIVE_EXPANSION[name]:
                if gg.name == "cx":
                    a, b = gg.qubits
                    v = v.reshape(2, 2)
                    if (a, b) == (0, 1):
                        v[1] = v[1][::-1]
                    else:
                        v[:, 0] = v[:, 0][::-1]
                    v = v.reshape(4)
                else:
                    m1 = np.array(gate_matrix(gg), dtype=complex)
                    w = gg.qubits[0]
                    v = v.reshape(2, 2)
                    v = m1 @ v if w == 0 else (m1 @ v.T).T
                    v = v.reshape(4)
            return v
        for col in range(4):
            e = np.zeros(4, dtype=complex)
            e[col] = 1.0
            M[:, col] = apply_col(e)
        return M
    raise ValueError(f"verify_large: unsupported 2q gate {name!r}")


def state_apply(circ, sv, deadline=None):
    """Apply circuit to a 2^n statevector (numpy).  `sv` is modified copy.
    `deadline` (absolute perf_counter time) bounds the application: a
    wide slow circuit declines instead of running for hours."""
    np = _np()
    from .linalg import gate_matrix
    from .mcx import expand_gate
    n = circ.num_qubits
    sv = sv.copy()
    for g in circ.ops:
        if deadline is not None and time.perf_counter() > deadline:
            raise ValueError("randomized verification exceeded its "
                             "wall-clock deadline")
        if len(g.qubits) == 1:
            m = np.asarray(gate_matrix(g), dtype=complex).reshape(2, 2)
            t = sv.reshape([2] * n)
            t = np.moveaxis(t, n - 1 - g.qubits[0], 0).reshape(2, -1)
            t = m @ t
            sv = np.moveaxis(t.reshape([2] * n), 0,
                             n - 1 - g.qubits[0]).ravel()
        elif len(g.qubits) == 2:
            M = _mat2q(g.name, g.params)
            a, b = g.qubits
            t = sv.reshape([2] * n)
            t = np.moveaxis(t, [n - 1 - a, n - 1 - b], [0, 1]).reshape(4, -1)
            t = M @ t
            sv = np.moveaxis(t.reshape([2] * n), [0, 1],
                             [n - 1 - a, n - 1 - b]).ravel()
        else:
            ex = expand_gate(g)
            if ex is None:
                raise ValueError(f"verify_large: unsupported gate {g.name!r}")
            sv = state_apply(type(circ)(circ.num_qubits, ex), sv, deadline)
    return sv


def states_agree(circ_a, circ_b, k: int = 32, tol: float = 1e-8, seed: int = 0,
                 deadline_s: float | None = None):
    """True iff both circuits map K random product states to identical
    outputs (within tol).  Probabilistically exact; K=32 gives failure
    probability far below hardware error rates.

    Refuses widths above ~30 qubits (statevector memory) and declines
    past `deadline_s` (wall clock) — the caller falls through, never
    gets a guessed answer."""
    np = _np()
    if circ_a.num_qubits != circ_b.num_qubits:
        return False
    n = circ_a.num_qubits
    if n > RAND_MAX_QUBITS:
        raise ValueError(
            f"randomized verification limited to {RAND_MAX_QUBITS} qubits "
            "(statevector memory); the caller declines beyond this")
    rng = np.random.default_rng(seed)
    deadline = time.perf_counter() + deadline_s if deadline_s else None
    for _ in range(k):
        psi = np.ones(1, dtype=complex)
        for _ in range(n):
            v = rng.normal(size=2) + 1j * rng.normal(size=2)
            v /= np.linalg.norm(v)
            psi = np.kron(psi, v)
        da = state_apply(circ_a, psi, deadline)
        db = state_apply(circ_b, psi, deadline)
        # phase-insensitive: circuits are equal up to global phase
        na = np.linalg.norm(da)
        nb = np.linalg.norm(db)
        ov = abs(np.vdot(da, db))
        if ov / (na * nb) < 1 - tol:
            return False
    return True


def optimize_large(circ, k: int = 32, tol: float = 1e-8, max_qubits: int = 30):
    """optimize_search without the dense prover, then the strongest
    available proof: exact algebraic (phase-polynomial), exact
    decision-diagram, and finally randomized K-state verification.

    Returns (circuit, status): "exact-proven (algebraic)" and
    "exact-proven (decision-diagram)" are machine-checked exact proofs;
    "exact" means K-state randomized verification passed;
    "rejected" means verification failed and the ORIGINAL circuit is
    returned instead (never ship unproven).
    """
    from . import optimize_search
    from .optimize import _UNVERIFIED
    if circ.num_qubits > max_qubits:
        raise ValueError(f"optimize_large supports <= {max_qubits} qubits")
    out = optimize_search(circ, verify=_UNVERIFIED)

    # exact tier 3: CX+diagonal fragment -> algebraic proof, any width
    from .phasepoly import phasepoly_equal
    try:
        pp = phasepoly_equal(circ, out)
    except Exception:
        pp = None
    if pp is True:
        return out, "exact-proven (algebraic)"
    if pp is False:
        return circ, "rejected"

    # exact tier 2: decision-diagram proof (node-budgeted AND
    # wall-clock-deadline-bounded; a blowup or a slow diagram is a
    # decline, never a wrong answer — same deadline contract as the
    # verify() cascade)
    from .dd import check_equivalent_dd, dd_fidelity, default_max_nodes
    import os as _os
    raw = _os.environ.get("COMPACTQ_DD_DEADLINE_S", "120")
    deadline = float(raw) if raw.strip() else None
    nodes = default_max_nodes()
    try:
        if check_equivalent_dd(circ, out, max_nodes=nodes,
                               deadline_s=deadline):
            return out, "exact-proven (decision-diagram)"
        if dd_fidelity(circ, out, max_nodes=nodes,
                       deadline_s=deadline) < 1.0 - 1e-7:
            return circ, "rejected"
    except Exception:
        pass  # budget / deadline exhausted / unsupported gate: fall to randomized

    try:
        raw2 = _os.environ.get("COMPACTQ_VERIFY_DEADLINE_S", "120")
        rand_deadline = float(raw2) if raw2.strip() else None
        if states_agree(circ, out, k=k, tol=tol, deadline_s=rand_deadline):
            return out, "exact"
        return circ, "rejected"
    except (MemoryError, ValueError):
        # statevector memory exhausted, width refused, or deadline hit:
        # could not decide — return the original (never ship unproven),
        # loudly
        return circ, "unverified (prover declined)"
