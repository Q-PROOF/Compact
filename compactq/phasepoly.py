"""Algebraic phase-polynomial prover — exact proof for the
CX + diagonal-rotation fragment at ANY qubit count.

A circuit over {cx, swap, rz, p, s, sdg, z, t, tdg, id} implements

    U |x>  =  e^{i phi(x)} . |A x>

where A is a GF(2) linear map (CNOT conjugation, SWAP permutation) and
phi(x) = sum over input-bit masks m of theta_m . [ <m, x> mod 2 = 1 ]
is a real linear combination of input parities (every diagonal rotation
adds its angle to the parity carried by its wire at that point; the
per-gate constants are global phases and drop out).

Two such circuits are equal up to global phase EXACTLY when their
linear maps and their parity-angle tables coincide — an algebraic
decision with exact GF(2) arithmetic; the only numerics are angle sums
compared at the house tolerance (1e-7, the same grade as the dense
prover).  This proves QFT / Trotter / phase-rotation circuits exactly
at widths where no dense matrix and no statevector can exist.

Gate conventions: rz(t) = diag(e^{-it/2}, e^{+it/2}) contributes +t to
its wire's parity (constant -t/2 is global); p(t) = diag(1, e^{+it})
contributes +t; s/sdg/z/t/tdg are p at +-pi/2, pi, +-pi/4.
"""
from __future__ import annotations

import math

from .circuit import Circuit

__all__ = ["phasepoly_key", "phasepoly_equal", "is_phasepoly_circuit",
           "PHASEPOLY_GATES"]

PHASEPOLY_GATES = frozenset(
    {"cx", "swap", "rz", "p", "s", "sdg", "z", "t", "tdg", "id"})

_FIXED_ANGLE = {
    "s": math.pi / 2,
    "sdg": -math.pi / 2,
    "z": math.pi,
    "t": math.pi / 4,
    "tdg": -math.pi / 4,
}


def is_phasepoly_circuit(circ: Circuit) -> bool:
    """True when every gate is in the CX + diagonal-rotation fragment."""
    return all(g.name in PHASEPOLY_GATES for g in circ.ops)


def phasepoly_key(circ: Circuit):
    """Return ((mask_per_wire...), {mask: angle}) for the circuit, or
    None when the circuit leaves the fragment.

    mask_per_wire[w] is the n-bit mask of input bits whose parity equals
    wire w's value; the phase table maps an input-parity mask to its
    accumulated rotation angle (global constants dropped).
    """
    n = circ.num_qubits
    masks = [1 << w for w in range(n)]
    phases: dict = {}
    for g in circ.ops:
        name = g.name
        if name == "cx":
            c, t = g.qubits
            masks[t] ^= masks[c]
        elif name == "swap":
            a, b = g.qubits
            masks[a], masks[b] = masks[b], masks[a]
        elif name == "id":
            continue
        elif name in ("rz", "p"):
            if not g.params:
                return None
            ang = float(g.params[0])
            if abs(ang) > 1e-13:
                m = masks[g.qubits[0]]
                phases[m] = phases.get(m, 0.0) + ang
        elif name in _FIXED_ANGLE:
            m = masks[g.qubits[0]]
            phases[m] = phases.get(m, 0.0) + _FIXED_ANGLE[name]
        else:
            return None
    return tuple(masks), phases


def phasepoly_equal(circ_a: Circuit, circ_b: Circuit,
                    tol: float = 1e-7) -> bool | None:
    """Decisive algebraic equivalence for fragment circuits.

    True / False when BOTH circuits are in the fragment (False is exact:
    differing GF(2) maps or parity tables are genuinely inequivalent);
    None when either circuit leaves the fragment — the caller falls
    through to the next prover in the cascade.
    """
    if circ_a.num_qubits != circ_b.num_qubits:
        return False
    ka = phasepoly_key(circ_a)
    if ka is None:
        return None
    kb = phasepoly_key(circ_b)
    if kb is None:
        return None
    masks_a, phases_a = ka
    masks_b, phases_b = kb
    if masks_a != masks_b:
        return False
    if set(phases_a) != set(phases_b):
        return False
    two_pi = 2.0 * math.pi
    for m, ang in phases_a.items():
        # phases are modulo 2*pi: resynthesis legitimately emits angles
        # shifted by whole turns, so compare the wrapped difference
        d = (ang - phases_b[m] + math.pi) % two_pi - math.pi
        if abs(d) > tol:
            return False
    return True
