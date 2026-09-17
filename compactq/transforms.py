"""Optimization passes.

Three cooperating passes, iterated to a fixpoint by ``optimize``:

1. ``peephole``         — folds every maximal run of 1-qubit gates on a wire
                          into at most 3 canonical gates (exact up to global
                          phase).  Cancels inverse pairs and merges rotations.
2. ``commute_cancel``   — slides provably-commuting gates out of the way, then
                          cancels self-inverse CX/CZ pairs it exposes.
3. ``swap_template``    — rewrites CX(a,b) CX(b,a) CX(a,b) into SWAP.

Every rewrite is *exact* (unitary preserved up to global phase); the
property tests prove this over thousands of random circuits.
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate
from .linalg import (gate_matrix, is_diagonal, mmul, mident, resynth,
                     same_up_to_phase)

# Gates we may fold into a named single gate after a run collapses.
_SINGLE = ("x", "y", "z", "h", "s", "sdg", "t", "tdg")


# --------------------------------------------------------------------- utils
def _folded_to_gates(prod, q: int):
    """Turn a folded 2x2 product into gate list on qubit q (exact up to phase)."""
    if prod is None:
        return None
    if mident(prod):
        return []
    for name in _SINGLE:
        from .linalg import NAMED_MATRICES
        if same_up_to_phase(prod, NAMED_MATRICES[name]):
            return [Gate(name, (), (q,))]
    gates = resynth(prod)
    if gates is None:
        return None
    return [Gate(g.name, g.params, (q,)) for g in gates]




# ------------------------------------------------- generalized diagonal slide
_DIAG_1Q = frozenset({"p", "rz", "s", "sdg", "z", "t", "tdg"})


def diag_slide(circ: Circuit) -> Circuit:
    """Slide diagonal 1-qubit gates through anything they commute with.

    A diagonal gate on wire q hops left over:
      * any gate not touching q (different wires always commute),
      * CX(c,t) where q is the control  (diagonals commute with the control),
      * CZ where q is either wire       (CZ is diagonal).

    This merges the control phases of adjacent CP blocks across block
    boundaries — the exact transformation that compresses QFT-style circuits.
    """
    ops = list(circ.ops)
    n = len(ops)
    if n < 2:
        return circ
    for _sweep in range(3):   # bounded: cross-wire diagonal pairs may re-swap
        moved = False
        i = n - 1
        while i > 0:
            g = ops[i]
            if len(g.qubits) != 1 or g.name not in _DIAG_1Q:
                i -= 1
                continue
            q = g.qubits[0]
            prev = ops[i - 1]
            hop = False
            if q not in prev.qubits:
                hop = True                       # different wires: commute
            elif prev.name == "cx" and prev.qubits[0] == q:
                hop = True                       # diagonal commutes with control
            elif prev.name == "cz":
                hop = True                       # diagonal commutes with CZ
            if hop:
                ops[i - 1], ops[i] = g, prev
                moved = True
            i -= 1
    return Circuit(circ.num_qubits, ops)


# ------------------------------------------------------------------ peephole
def peephole(circ: Circuit) -> Circuit:
    """Fold runs of 1-qubit gates per wire into <= 3 gates.

    A folded run replaces the original gates only when it is *strictly
    shorter* — the optimizer never grows a circuit through resynthesis.
    """
    out = []
    runs: dict = {}  # qubit -> [product matrix, op count, original ops]

    def flush(q: int) -> None:
        if q not in runs:
            return
        prod, count, orig = runs.pop(q)
        gates = _folded_to_gates(prod, q)
        if gates is not None and len(gates) < count:
            out.extend(gates)
        else:
            out.extend(orig)

    for g in circ.ops:
        if len(g.qubits) == 1:
            q = g.qubits[0]
            m = gate_matrix(g)
            if q in runs:
                # circuit semantics are left-multiplication: U = G_k ... G_1,
                # so each newly-seen gate goes on the LEFT of the product.
                runs[q][0] = mmul(m, runs[q][0])
                runs[q][1] += 1
                runs[q][2].append(g)
            else:
                runs[q] = [m, 1, [g]]
        else:
            for q in g.qubits:
                flush(q)
            out.append(g)
    for q in list(runs):
        flush(q)
    return Circuit(circ.num_qubits, out)


# ------------------------------------------------------- commute and cancel
def _cx_rewrite(pair, accum):
    """Rebuild the gate sequence between a cancelled CX pair."""
    c, t = pair
    out = []
    for h in accum:
        if len(h.qubits) != 1:
            out.append(h)
            continue
        q = h.qubits[0]
        if q == t and (h.name == "x" or h.name == "rx"):
            out.append(h)  # CX X_t CX = X_t  /  CX RX_t CX = RX_t
        elif q == c and is_diagonal(h):
            out.append(h)  # commutes with control
        elif q == c and h.name == "x":
            out.append(Gate("x", (), (c,)))
            out.append(Gate("x", (), (t,)))  # CX X_c CX = X_c X_t
        else:
            out.append(h)  # unrelated wire
    return out


def _cz_rewrite(pair, accum):
    c, t = pair
    out = []
    for h in accum:
        if len(h.qubits) != 1:
            out.append(h)
            continue
        q = h.qubits[0]
        if is_diagonal(h):
            out.append(h)  # diagonal commutes with CZ on either wire
        elif h.name == "x" and q == c:
            out.append(Gate("x", (), (c,)))
            out.append(Gate("z", (), (t,)))  # CZ X_c CZ = X_c Z_t
        elif h.name == "x" and q == t:
            out.append(Gate("z", (), (c,)))
            out.append(Gate("x", (), (t,)))  # CZ X_t CZ = Z_c X_t
        else:
            out.append(h)
    return out


def commute_cancel(circ: Circuit) -> Circuit:
    """Cancel self-inverse CX/CZ pairs across provably-commuting gates."""
    ops = circ.ops
    changed = True
    while changed:
        changed = False
        for i in range(len(ops)):
            g = ops[i]
            if g.name not in ("cx", "cz"):
                continue
            c, t = g.qubits
            accum = []
            j = i + 1
            matched = False
            while j < len(ops):
                h = ops[j]
                if len(h.qubits) == 1:
                    q = h.qubits[0]
                    if q == c and g.name == "cx":
                        if is_diagonal(h) or h.name == "x":
                            accum.append(h); j += 1; continue
                        break
                    if q == t and g.name == "cx":
                        if h.name in ("x", "rx"):
                            accum.append(h); j += 1; continue
                        break
                    if g.name == "cz":
                        if h.name == "x" or is_diagonal(h):
                            accum.append(h); j += 1; continue
                        break
                    accum.append(h); j += 1; continue  # unrelated wire
                hq = set(h.qubits)
                if len(hq) == 2:
                    if h.name == "cz" and hq == {c, t} and g.name == "cz":
                        matched = True
                        break
                    # CX must match in the SAME direction; a reversed CX pair
                    # composes to SWAP, not identity.
                    if h.name == "cx" and g.name == "cx" and h.qubits == g.qubits:
                        matched = True
                        break
                break
            if matched:
                rewrite = _cx_rewrite((c, t), accum) if g.name == "cx" else _cz_rewrite((c, t), accum)
                ops = ops[:i] + rewrite + ops[j + 1:]
                changed = True
                break
    return Circuit(circ.num_qubits, ops)


# ----------------------------------------------------------- 1q phase slide
def slide_1q(circ: Circuit) -> Circuit:
    """Slide 1-qubit gates left across 2-qubit gates they provably commute with.

    * Diagonal gates on a CX/CZ *control* wire commute with the 2-qubit gate.
    * Diagonal gates on either wire of a CZ commute with it.
    * X-type gates on a CX *target* wire commute with it.

    This pulls redundant phases out of 2-qubit blocks so `peephole` can
    absorb them into neighbouring runs.
    """
    ops = list(circ.ops)
    changed = True
    while changed:
        changed = False
        for i in range(len(ops) - 1):
            g, h = ops[i], ops[i + 1]
            if len(h.qubits) != 1:
                continue
            if g.name == "cx":
                c, t = g.qubits
                if (h.qubits[0] == c and is_diagonal(h)) or \
                   (h.qubits[0] == t and h.name in ("x", "rx")):
                    ops[i], ops[i + 1] = h, g
                    changed = True
                    break
            elif g.name == "cz" and is_diagonal(h) and h.qubits[0] in g.qubits:
                ops[i], ops[i + 1] = h, g
                changed = True
                break
    return Circuit(circ.num_qubits, ops)


# --------------------------------------------------------------- swap rule
def swap_template(circ: Circuit) -> Circuit:
    """CX(a,b) CX(b,a) CX(a,b)  ->  SWAP(a, b)."""
    ops = circ.ops
    out = []
    i = 0
    while i < len(ops):
        if (i + 2 < len(ops)
                and ops[i].name == ops[i + 1].name == ops[i + 2].name == "cx"
                and ops[i].qubits[0] == ops[i + 1].qubits[1]
                and ops[i].qubits[1] == ops[i + 1].qubits[0]
                and ops[i + 2].qubits == ops[i].qubits):
            a, b = ops[i].qubits
            out.append(Gate("swap", (), (a, b)))
            i += 3
        else:
            out.append(ops[i])
            i += 1
    return Circuit(circ.num_qubits, out)
