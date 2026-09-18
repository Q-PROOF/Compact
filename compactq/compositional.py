"""T4 compositional verification for block-structured circuits.

Soundness argument: if the active qubits of a circuit partition into
disjoint components (every gate acts inside one component), the circuit's
unitary factorizes as the tensor product of per-component unitaries.  Two
such circuits are equivalent (up to global phase) iff every matched
component pair is equivalent — so each block can be dense-verified
independently at its OWN width, and no 2^total_qubits unitary is ever
built.  This extends the equivalence-proven regime to wide
block-structured circuits (e.g. a 100-qubit circuit of 4-qubit blocks).
"""
from __future__ import annotations

from .circuit import Circuit, Gate

__all__ = ["active_components", "verify_compositional"]


def active_components(circ: Circuit):
    """Split `circ` into per-component subcircuits over its ACTIVE qubits.

    Returns a list of (frozenset_of_qubits, subcircuit) sorted by qubits.
    Idle qubits carry identity and are ignored (both sides must agree on
    the total qubit count, which verify() already checks).
    """
    ops = circ.ops
    if not ops:
        return []
    parent = {q: q for g in ops for q in g.qubits}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in ops:
        r0 = find(g.qubits[0])
        for q in g.qubits[1:]:
            parent[find(q)] = r0

    groups: dict[int, list] = {}
    for g in ops:
        groups.setdefault(find(g.qubits[0]), []).append(g)

    comps = []
    for _root, gops in groups.items():
        qs = sorted({q for g in gops for q in g.qubits})
        remap = {q: i for i, q in enumerate(qs)}
        sub = Circuit(len(qs),
                      [Gate(g.name, g.params,
                            tuple(remap[q] for q in g.qubits)) for g in gops])
        comps.append((frozenset(qs), sub))
    comps.sort(key=lambda t: sorted(t[0]))
    return comps


def verify_compositional(a: Circuit, b: Circuit, tol: float = 1e-7):
    """Compositional (tier 4) equivalence check for block-structured pairs.

    Returns None when the method cannot decide (single component, mismatched
    component structure, or a block wider than the dense prover's limit) —
    callers fall back to a weaker tier.  Otherwise returns an evidence dict;
    `equivalent: False` names the failing block's qubits.
    """
    from .equivalence import check_equivalent, _MAX_QUBITS

    ca = active_components(a)
    cb = active_components(b)
    if len(ca) < 2 or len(cb) < 2:
        return None
    if [qs for qs, _ in ca] != [qs for qs, _ in cb]:
        return None  # different block structure: out of scope for T4
    for (qs, sub_a), (_qs2, sub_b) in zip(ca, cb):
        if len(qs) > _MAX_QUBITS:
            return None  # block itself beyond the dense prover
        if not check_equivalent(sub_a, sub_b, tol):
            return {"equivalent": False, "tier": 4,
                    "method": "compositional_disjoint_blocks",
                    "blocks": len(ca),
                    "failing_block_qubits": sorted(qs),
                    "global_phase_ignored": True}
    return {"equivalent": True, "tier": 4,
            "method": "compositional_disjoint_blocks",
            "blocks": len(ca), "global_phase_ignored": True}
