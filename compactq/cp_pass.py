"""Commutative cancellation of CX pairs around target phases, producing
controlled-phase (cp) gates.

Exact identity (verified numerically; equal up to global phase e^{i t/2}):
  CX(c,t) . RZ_t(t) . CX(c,t)  =  RZ_c(t) . RZ_t(t) . CP(-2t)

since the CX pair conjugates a target rz into the ZZ rotation
RZZ(t) = exp(-i(t/2) Zc Zt).  So [cx, rz_t(θ), cx]  ->  [rz_c(θ), rz_t(θ),
CP(-2θ)] — two CX become one CP.  Adjacent CPs on the same pair then merge
(CP(θ)·CP(φ) = CP(θ+φ)), and CP(0) drops out.  For QFT/QAOA/PEA-style
circuits this halves — or better — the two-qubit count.
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate
from .symbolic import Sym, is_zero as _is_zero

_DIAG = frozenset({"p", "rz", "s", "sdg", "z", "t", "tdg"})


def _is_phase_diag(g: Gate) -> bool:
    """True if g is a diagonal 1-qubit gate expressible as RZ(angle)."""
    return g.name in _DIAG


def _diag_angle(g: Gate) -> float:
    return g.params[0] if g.name in ("rz", "p") else _PHASE_OF[g.name]


_PHASE_OF = {"s": math.pi / 2, "sdg": -math.pi / 2, "z": math.pi,
             "t": math.pi / 4, "tdg": -math.pi / 4}


def _wrap(x):
    if isinstance(x, Sym):
        return x  # symbolic angles stay unwrapped (exact for RZ/CP)
    return (x + math.pi) % (2 * math.pi) - math.pi


def cancel_cx_through_diagonal(circ: Circuit, max_span: int = 64) -> Circuit:
    """Find matched CX pairs [cx(c,t), ..., cx(c,t)] and rewrite them via
      CX(c,t) . RZ_t(Θ) . CX(c,t)  =  RZ_c(Θ) . RZ_t(Θ) . CP(-2Θ)
    (up to global phase).  Between the pair, every gate must be either
      - a 1-qubit diagonal on t (angles fold into Θ) or on c (they commute
        through the CX and fold into the RZ_c),
      - a 2-qubit diagonal on exactly {c,t} (cp/cz, angles fold into Φ), or
      - a gate on wires disjoint from {c,t} (commutes with the pair; kept
        in front of the rewrite in original order).
    Any other gate on {c,t} aborts the scan for this pair.

    With Φ the sum of the cp/cz phases, the full rewrite is
      rz_c(Σa + Θ + Φ/2) . rz_t(Θ - Φ/2) . cp(-2Θ - Φ)
    and when all three angles vanish the pair disappears entirely.
    Replaces two CX with one CP.  This is what compresses QAOA/QFT/PEA
    phase ladders and relative-phase adder (T-sandwiched CX) patterns.
    """
    ops = list(circ.ops)
    out = []
    i = 0
    n = len(ops)
    while i < n:
        g = ops[i]
        if g.name == "cx" and i + 2 < n:
            c, t = g.qubits
            theta_t = 0.0
            theta_c = 0.0
            phi_cp = 0.0
            keep = []
            j = i + 1
            matched = False
            while j < n and j - i <= max_span:
                m = ops[j]
                if m.name == "cx" and m.qubits == (c, t):
                    matched = True
                    break
                if len(m.qubits) == 2 and m.name in ("cp", "cz") \
                        and set(m.qubits) == {c, t}:
                    phi_cp += m.params[0] if m.name == "cp" else math.pi
                elif len(m.qubits) == 1 and m.qubits[0] in (c, t) \
                        and m.name in _DIAG:
                    # rz/p with no params (rz(0)) are legal and fold to 0
                    if m.name in ("rz", "p"):
                        ang = m.params[0] if m.params else 0.0
                    else:
                        ang = _PHASE_OF[m.name]
                    if m.qubits[0] == t:
                        theta_t += ang
                    else:
                        theta_c += ang
                elif set(m.qubits).isdisjoint((c, t)):
                    keep.append(m)
                else:
                    break
                j += 1
            if matched:
                theta_t_w = _wrap(theta_t)
                phi_w = _wrap(phi_cp)
                # CX.CP(phi).CX = RZ_c(phi).CP(-phi) (up to phase), so the
                # folded 2q diagonals add Phi to the control rz and -Phi to
                # the cp; the target rz is untouched
                rz_c = _wrap(theta_c + theta_t_w + phi_w)
                rz_t = _wrap(theta_t_w)
                cp = _wrap(-2.0 * theta_t_w - phi_w)
                if _is_zero(rz_c) and _is_zero(rz_t) and _is_zero(cp):
                    # everything folds away: the pair vanishes
                    out.extend(keep)
                    i = j + 1
                    continue
                out.extend(keep)
                if not _is_zero(rz_c):
                    out.append(Gate("rz", (rz_c,), (c,)))
                if not _is_zero(rz_t):
                    out.append(Gate("rz", (rz_t,), (t,)))
                if not _is_zero(cp):
                    out.append(Gate("cp", (cp,), (c, t)))
                i = j + 1
                continue
        out.append(g)
        i += 1
    return Circuit(circ.num_qubits, out)


def merge_cp(circ: Circuit) -> Circuit:
    """Merge adjacent/same-pair CP gates: CP(θ)·CP(φ) = CP(θ+φ); drop CP(0)."""
    ops = list(circ.ops)
    # merge consecutive cp on the same pair
    out = []
    i = 0
    n = len(ops)
    while i < n:
        g = ops[i]
        if g.name == "cp":
            theta = g.params[0]
            j = i + 1
            while j < n and ops[j].name == "cp" and ops[j].qubits == g.qubits:
                theta += ops[j].params[0]
                j += 1
            theta = _wrap(theta)
            if not _is_zero(theta):
                out.append(Gate("cp", (theta,), g.qubits))
            i = j
        else:
            out.append(g)
            i += 1
    return Circuit(circ.num_qubits, out)
