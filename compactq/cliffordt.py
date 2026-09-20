"""Clifford+T normal-form optimization.

For circuits whose gates are all in the Clifford+T set (H, X, S-family,
CX, SWAP, CZ, T-family, RZ/P at multiples of pi/4, CP at multiples of
pi/4), the circuit factors into an alternating sequence

    C0 . P0 . C1 . P1 . ... . Ck

of Clifford layers C_i (tableau-exact) and diagonal phase layers P_i
(parity polynomials over {0,1}^n with angles in pi/4 units).  Each P_i is
re-synthesized as a shared parity network (CNOT ladders reused across
wires) and each C_i as an AG tableau circuit *with free output wire
permutations* - a permutation chosen for one block is absorbed by
relabeling everything downstream, so SWAPs are only paid once, at the
very end of the circuit.

Exactness: the layer partition preserves gate order, every Clifford
replacement carries a tableau proof, the parity networks are exact, and
the final circuit is re-verified against the input (dense prover at <= 8
qubits) before it can replace anything.  Circuits containing any gate
outside the Clifford+T set return None (caller keeps the current
pipeline).
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate
from .stabilizer import Tableau
from .parity import _parity_network
from .clifford import _synth_gates_from_tableau

_Q4 = math.pi / 4
_TOL = 1e-9


def _classify(g: Gate):
    """Classify a gate into ("C", gate) / ("P", [phase terms]) / None.

    Phase terms are (wire_mask, rz_angle) contributions.
    """
    n = len(g.qubits)
    name = g.name
    th = g.params[0] if g.params else 0.0
    if n == 1:
        w = g.qubits[0]
        if name in ("h", "x", "s", "sdg", "z"):
            return "C", [g]
        if name == "id":
            return "C", []
        if name in ("t",):
            return "P", [(1 << w, _Q4)]
        if name == "tdg":
            return "P", [(1 << w, -_Q4)]
        if name in ("rz", "p"):
            k = int(round(th / _Q4))
            if abs(th - k * _Q4) > _TOL:
                return None, None
            r = k % 4
            if r == 0:
                return "C", []
            if r == 2:
                # rz(pi/2) or rz(-pi/2): Clifford (S / Sdg up to phase)
                out = "s" if r == 2 and k > 0 else "sdg"
                return "C", [Gate(out, (), (w,))]
            ang = _Q4 if r == 1 else -_Q4
            return "P", [(1 << w, ang)]
        return None, None
    if n == 2:
        a, b = g.qubits
        if name == "cx":
            return "C", [g]
        if name == "swap":
            return "C", [Gate("cx", (), (a, b)), Gate("cx", (), (b, a)),
                         Gate("cx", (), (a, b))]
        if name == "cz":
            # CZ = (H x H) . CX . (H x H) up to nothing (diagonal conj trick)
            return "C", [Gate("h", (), (a,)), Gate("h", (), (b,)),
                         Gate("cx", (), (a, b)),
                         Gate("h", (), (a,)), Gate("h", (), (b,))]
        if name == "cp":
            k = int(round(th / _Q4))
            if abs(th - k * _Q4) > _TOL:
                return None, None
            r = k % 4
            if r == 0:
                return "C", []
            if r == 2:
                # CP(pi) = CZ exactly
                return "C", [Gate("h", (), (a,)), Gate("h", (), (b,)),
                             Gate("cx", (), (a, b)),
                             Gate("h", (), (a,)), Gate("h", (), (b,))]
            # CP(th) phase polynomial (exact): {a: -th/2, b: -th/2, ab: +th}
            return "P", [((1 << a), -th / 2), ((1 << b), -th / 2),
                         ((1 << a) | (1 << b), th)]
        return None, None
    return None, None


def _segments(circ: Circuit):
    """Yield alternating ("C", gates) / ("P", terms) segments, or None."""
    segs = []
    cur_kind = None
    cur_gates = []
    cur_terms = {}
    for g in circ.ops:
        kind, payload = _classify(g)
        if kind is None:
            return None
        if kind == "C":
            if cur_kind == "P":
                segs.append(("P", dict(cur_terms)))
                cur_terms = {}
                cur_kind = None
            cur_kind = "C"
            cur_gates.extend(payload)
        else:
            if cur_kind == "C":
                segs.append(("C", list(cur_gates)))
                cur_gates = []
                cur_kind = None
            cur_kind = "P"
            for mask, ang in payload:
                cur_terms[mask] = cur_terms.get(mask, 0.0) + ang
    if cur_kind == "C":
        segs.append(("C", cur_gates))
    elif cur_kind == "P":
        segs.append(("P", cur_terms))
    return segs


def _tableau_of_gates(gates, n):
    circ = Circuit(n, gates)
    tab = Tableau(n)
    for g in circ.ops:
        if len(g.qubits) == 1:
            from .clifford import _is_clifford_1q
            table = _is_clifford_1q(g)
            if table is None:
                return None
            tab.apply_1q_table(g.qubits[0], table)
        elif g.name == "cx":
            tab.apply_cx(g.qubits[0], g.qubits[1])
        else:
            return None
    return tab


def _synth_count(gates, n):
    return sum(1 for g in gates if len(g.qubits) == 2)


def cliffordt_pass(circ: Circuit, verify: bool = True,
                   absorb_perms: bool = True):
    """Clifford+T normal form; returns an exact Circuit or None.

    Accepts only circuits over the Clifford+T gate set (see module doc).
    The rebuilt circuit is verified against the input (dense prover) and
    returned only when strictly smaller on 2-qubit count.
    """
    segs = _segments(circ)
    if segs is None or len(segs) == 0:
        return None
    n = circ.num_qubits
    perm = list(range(n))            # original wire -> current wire
    out = []
    total_cx = 0

    for kind, payload in segs:
        if kind == "C":
            if not payload:
                continue
            # relabel gates into current wire positions
            rel = [Gate(g.name, g.params, tuple(perm[w] for w in g.qubits))
                   for g in payload]
            tab = _tableau_of_gates(rel, n)
            if tab is None:
                return None
            base = _synth_gates_from_tableau(tab)
            if base is None:
                return None
            best_gates, best_cx, best_pi = base, _synth_count(base, n), None
            wires = sorted({w for g in rel for w in g.qubits})
            if absorb_perms and len(wires) <= 5 and best_cx > 0:
                from itertools import permutations
                ident = list(range(n))
                for pi_part in permutations(wires):
                    if list(pi_part) == wires:
                        continue
                    mapping = dict(zip(wires, pi_part))
                    rel2 = [Gate(g.name, g.params,
                                 tuple(mapping[w] for w in g.qubits))
                            for g in rel]
                    tab2 = _tableau_of_gates(rel2, n)
                    if tab2 is None:
                        continue
                    g2 = _synth_gates_from_tableau(tab2)
                    if g2 is None:
                        continue
                    c2 = _synth_count(g2, n)
                    if c2 < best_cx:
                        best_gates, best_cx, best_pi = g2, c2, mapping
            out.extend(best_gates)
            total_cx += best_cx
            if best_pi is not None:
                # emitted gates = P_pi . G_block . P_pi^dag: downstream relabels by pi
                perm = [best_pi[perm[w]] if perm[w] in best_pi else perm[w]
                        for w in range(n)]
        else:
            if not payload:
                continue
            # relabel term masks into current wire positions
            terms = {}
            for mask, ang in payload.items():
                new_mask = 0
                for w in range(n):
                    if (mask >> w) & 1:
                        new_mask |= 1 << perm[w]
                terms[new_mask] = terms.get(new_mask, 0.0) + ang
            gates = _parity_network(terms, n)
            if gates is None:
                return None
            out.extend(gates)
            total_cx += sum(1 for g in gates if len(g.qubits) == 2)

    # fix the residual wire permutation with explicit SWAPs
    swaps = _permutation_swaps(perm)
    for a, b in swaps:
        out.extend([Gate("cx", (), (a, b)), Gate("cx", (), (b, a)),
                    Gate("cx", (), (a, b))])
        total_cx += 3

    if total_cx >= circ.two_qubit_count():
        return None
    rebuilt = Circuit(n, out)
    if verify:
        from .equivalence import check_equivalent
        if not check_equivalent(circ, rebuilt):
            return None
    return rebuilt


def _permutation_swaps(perm):
    """SWAP transpositions restoring perm to identity (perm[w] = current
    wire holding original wire w's content)."""
    pos = list(perm)
    swaps = []
    for w in range(len(pos)):
        while pos[w] != w:
            cur = pos[w]
            # swap contents of wires w and cur
            pos[w], pos[cur] = pos[cur], pos[w]
            swaps.append((w, cur))
    return swaps
