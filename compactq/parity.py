"""Phase-polynomial re-synthesis of diagonal cores.

A maximal window of CX gates interleaved with RZ-family diagonal gates
implements a phase polynomial  exp(i . sum_k theta_k . f_k(x))  where the
f_k are linear parity functions of the basis bits.  This module extracts
the parity terms from such a window (simulating CX parity propagation;
CP gates are first expanded into cx+rz exactly) and re-synthesizes them
with a parity network found by BFS over GF(2) parity states plus a full
uncompute stage, so the CX product cancels exactly and only the diagonal
phase remains.  Every replacement is verified against the window's own
unitary (<= 6 qubits) before it is accepted.
"""
from __future__ import annotations

from collections import deque

from .circuit import Circuit, Gate
from .equivalence import check_equivalent
from .linalg import _PHASE_ANGLE

_DIAG = ("rz", "p", "s", "sdg", "z", "t", "tdg")


def _angle_of(g: Gate) -> float:
    if g.name in ("rz", "p"):
        return g.params[0] if g.params else 0.0
    fn = _PHASE_ANGLE[g.name]
    return fn(g.params) if callable(fn) else fn


def collect_parity_windows(circ: Circuit, min_cx: int = 2):
    """Maximal contiguous windows composed only of CX and diagonal 1q gates.

    Returns (start, end, terms, n_qubits) where terms maps a bitmask support
    (bit w = wire w) to a total angle.  CP gates are expanded exactly into
    cx+rz before extraction.
    """
    ops = []
    for g in circ.ops:
        if g.name == "cp":
            a, b = g.qubits
            th = g.params[0]
            ops.extend([
                Gate("rz", (-th / 2,), (a,)),
                Gate("rz", (-th / 2,), (b,)),
                Gate("cx", (), (a, b)),
                Gate("rz", (th,), (b,)),
                Gate("cx", (), (a, b)),
            ])
        else:
            ops.append(g)

    windows = []
    i = 0
    n_ops = len(ops)
    n_q = circ.num_qubits
    while i < n_ops:
        g = ops[i]
        is_diag = len(g.qubits) == 1 and g.name in _DIAG
        is_cx = g.name == "cx"
        if not (is_diag or is_cx):
            i += 1
            continue
        j = i
        parity = [1 << w for w in range(n_q)]
        terms = {}
        n_cx = 0
        while j < n_ops:
            g = ops[j]
            if g.name == "cx":
                c, t = g.qubits
                parity[t] ^= parity[c]
                n_cx += 1
                j += 1
            elif len(g.qubits) == 1 and g.name in _DIAG:
                w = g.qubits[0]
                terms[parity[w]] = terms.get(parity[w], 0.0) + _angle_of(g)
                j += 1
            else:
                break
        if n_cx >= min_cx and terms:
            windows.append((i, j, terms, n_q))
            i = j
        else:
            i = j if j > i else i + 1
    return windows


def _covered(L, supports):
    ls = set(L)
    return all(s in ls for s, _ in supports)


def _parity_network(terms: dict, n: int, max_states: int = 4_000):
    """BFS a CX sequence making every support appear as some wire parity.

    Returns the circuit [CX path with interleaved rz] + [uncompute CXs],
    or None when the search budget is exhausted.  Uses the Rust kernel
    when available; the Python implementation stays as the reference and
    fallback.
    """
    try:
        import compactq_native as _nat
        if hasattr(_nat, "parity_network"):
            tl = [(s, a) for s, a in terms.items()]
            ops = _nat.parity_network(n, tl, max_states)
            out: list = []
            for kind, a, b in ops:
                if kind == 0:
                    out.append(Gate("rz", (b,), (a,)))
                else:
                    out.append(Gate("cx", (), (a, int(b))))
            return out
    except Exception:
        pass
    supports = [(s, a) for s, a in terms.items() if s and abs(a) > 1e-12]
    if not supports:
        return None
    start = tuple(1 << w for w in range(n))

    if _covered(start, supports):
        out = []
        for s, a in supports:
            w = next(w for w in range(n) if start[w] == s)
            out.append(Gate("rz", (a,), (w,)))
        return out

    seen = {start}
    frontier = deque([start])
    parent = {}
    goal = None
    while frontier:
        state = frontier.popleft()
        if state != start and _covered(state, supports):
            goal = state
            break
        if len(seen) >= max_states:
            break
        for c in range(n):
            for t in range(n):
                if c == t:
                    continue
                nl = list(state)
                nl[t] ^= nl[c]
                ns = tuple(nl)
                if ns in seen:
                    continue
                seen.add(ns)
                parent[ns] = (state, c, t)
                frontier.append(ns)
                if _covered(ns, supports):
                    goal = ns
                    frontier.clear()
                    break
    if goal is None:
        return None

    path = []
    cur = goal
    while cur != start:
        prev, c, t = parent[cur]
        path.append((c, t))
        cur = prev
    path.reverse()

    # replay the path, emitting rz as soon as its parity becomes available;
    # the uncompute reverses only the CX prefix actually emitted
    L = list(start)
    out = []
    cx_emitted = []
    pending = dict(supports)
    pi = 0
    while pending:
        progressed = False
        for s in [s for s in pending]:
            w = next((w for w in range(n) if L[w] == s), None)
            if w is not None:
                out.append(Gate("rz", (pending.pop(s),), (w,)))
                progressed = True
        if not pending:
            break
        if pi < len(path):
            c, t = path[pi]
            out.append(Gate("cx", (), (c, t)))
            cx_emitted.append((c, t))
            L[t] ^= L[c]
            pi += 1
            progressed = True
        if not progressed:
            return None
    for c, t in reversed(cx_emitted):
        out.append(Gate("cx", (), (c, t)))
    return out


def parity_pass(circ: Circuit) -> Circuit:
    """Re-synthesize diagonal cores with parity networks (exact, verified).

    Windows are only replaced when the circuit has <= 6 qubits so the
    replacement can be proven with a full-unitary fidelity check.
    """
    if circ.num_qubits > 6:
        return circ
    windows = collect_parity_windows(circ)
    if not windows:
        return circ
    ops = circ.ops
    out = []
    prev = 0
    changed = False
    for (start, end, terms, _n) in windows:
        out.extend(ops[prev:start])
        prev = end
        block_circ = Circuit(circ.num_qubits, ops[start:end])
        synth = _parity_network(dict(terms), circ.num_qubits)
        if synth is None:
            out.extend(ops[start:end])
            continue
        cand = Circuit(circ.num_qubits, synth)
        if not check_equivalent(block_circ, cand):
            out.extend(ops[start:end])
            continue
        if ((cand.two_qubit_count(), len(cand), cand.depth())
                < (block_circ.two_qubit_count(), len(block_circ),
                   block_circ.depth())):
            out.extend(synth)
            changed = True
        else:
            out.extend(ops[start:end])
    out.extend(ops[prev:])
    if not changed:
        return circ
    return Circuit(circ.num_qubits, out)

