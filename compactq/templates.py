"""Verified reordering templates for gate-level rewrites.

These templates REORDER gates (never add or remove) to expose
cancellation/merge opportunities that the standard passes cannot see
because a non-commuting-looking gate sits between two commuting ones:

  T3  cx(a,c) cx(b,c)   <->  cx(b,c) cx(a,c)    (same target: commute)
  T4  cx(a,b) cx(a,c)   <->  cx(a,c) cx(a,b)    (same control: commute)
  T5  diag1q(w) cz(a,b) <->  cz(a,b) diag1q(w)  (both diagonal: commute)
  T6  cx(a,b) z(b)      ->  z(a) z(b) cx(a,b)    (Z_target conjugates to
                                                   Z_control x Z_target)

T3/T4/T5 are exact reorderings; swaps only move toward a canonical order
so the fixpoint loop terminates.  T6 adds one gate but moves Z(b) left of
the CX where it can merge with adjacent phase gates; peephole makes the
net count neutral when a merge fires.

Every template is asserted against its matrix identity in the test suite,
and `template_pass` re-verifies the rewritten circuit against the original
(<= 6 qubits) before accepting it.
"""
from __future__ import annotations

from .circuit import Circuit, Gate
from .equivalence import check_equivalent


def _is(g: Gate, name: str):
    return g.name == name


def _rewrite(circ: Circuit) -> Circuit:
    ops = list(circ.ops)
    changed_any = False
    changed = True
    while changed:
        changed = False
        i = 0
        n = len(ops)
        while i < n:
            g = ops[i]
            if _is(g, "cx") and i + 1 < n and _is(ops[i + 1], "cx"):
                # T3: same target, different controls - canonical order:
                # smaller control first
                if g.qubits[1] == ops[i + 1].qubits[1] \
                        and g.qubits[0] > ops[i + 1].qubits[0]:
                    ops[i], ops[i + 1] = ops[i + 1], ops[i]
                    changed = changed_any = True
                    i += 1
                    continue
                # T4: same control, different targets - canonical order:
                # smaller target first
                if g.qubits[0] == ops[i + 1].qubits[0] \
                        and g.qubits[1] > ops[i + 1].qubits[1]:
                    ops[i], ops[i + 1] = ops[i + 1], ops[i]
                    changed = changed_any = True
                    i += 1
                    continue
            # T6: cx(a,b) followed by z(b)
            if _is(g, "cx") and i + 1 < n and _is(ops[i + 1], "z") \
                    and ops[i + 1].qubits[0] == g.qubits[1]:
                c, t = g.qubits
                ops = ops[:i] + [Gate("z", (), (c,)), Gate("z", (), (t,)),
                                 Gate("cx", (), (c, t))] + ops[i + 2:]
                changed = changed_any = True
                n = len(ops)
                i += 2
                continue
            # T5: diagonal 1q through CZ on either wire
            if _is(g, "cz") and i + 1 < n and len(ops[i + 1].qubits) == 1 \
                    and ops[i + 1].name in ("rz", "p", "s", "sdg", "z", "t", "tdg") \
                    and ops[i + 1].qubits[0] in g.qubits:
                ops[i], ops[i + 1] = ops[i + 1], ops[i]
                changed = changed_any = True
                i += 1
                continue
            i += 1
    if not changed_any:
        return circ
    return Circuit(circ.num_qubits, ops)


def template_pass(circ: Circuit) -> Circuit:
    """Apply reordering templates until fixpoint; verify for <= 6 qubits."""
    verify = circ.num_qubits <= 6
    out = _rewrite(circ)
    if out is circ:
        return circ
    if verify and not check_equivalent(circ, out):
        return circ
    return out
