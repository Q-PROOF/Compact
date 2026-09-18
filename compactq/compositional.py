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

__all__ = ["active_components", "verify_compositional", "verify_segmented"]


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


def verify_segmented(a: Circuit, b: Circuit, cuts_a, cuts_b,
                     tol: float = 1e-7):
    """T4.1 sequential-segment compositional proof.

    Soundness: for ANY valid consecutive segmentations of `a` and `b`, if
    every segment pair (A_i, B_i) is equivalent up to global phase then
    U_a = A_m···A_1 and U_b = B_m···B_1 are equivalent, because the per-
    segment phase factors commute with everything.  Each segment pair is
    dense-verified on its ACTIVE qubits only — so the method scales when
    the caller can name narrow segments (e.g. known Trotter layers,
    barrier boundaries, or a compiler's rewrite correspondence).

    IMPORTANT: the verdict is one-directional.  All segments matching ⇒
    `equivalent: True` (sound).  Any segment pair mismatching is only
    evidence that THIS segmentation does not align the two circuits — it
    is NOT a proof of global inequivalence, so the method returns
    `equivalent: None` (undecided) instead of False.

    `cuts_a` / `cuts_b` are strictly increasing gate-index cut positions
    (0 < cut < len(ops)).  Returns None when a segment pair's active
    width exceeds the dense prover, when the active qubit sets of a pair
    differ, or when the cut lists are malformed — callers fall back.
    """
    from .equivalence import check_equivalent, _MAX_QUBITS

    cuts_a = sorted(int(c) for c in cuts_a)
    cuts_b = sorted(int(c) for c in cuts_b)
    if len(cuts_a) != len(cuts_b):
        return None
    if (any(not 0 < c < len(a.ops) for c in cuts_a)
            or any(not 0 < c < len(b.ops) for c in cuts_b)):
        return None

    def segments(ops, cuts):
        bounds = [0] + cuts + [len(ops)]
        return [ops[bounds[i]:bounds[i + 1]] for i in range(len(bounds) - 1)]

    segs_a = segments(a.ops, cuts_a)
    segs_b = segments(b.ops, cuts_b)
    for i, (sa_ops, sb_ops) in enumerate(zip(segs_a, segs_b), start=1):
        qa = sorted({q for g in sa_ops for q in g.qubits})
        qb = sorted({q for g in sb_ops for q in g.qubits})
        if qa != qb:
            return None  # segments act on different qubits: not composable
        if not qa:
            continue  # both empty
        if len(qa) > _MAX_QUBITS:
            return None  # segment too wide for the dense prover
        remap = {q: j for j, q in enumerate(qa)}

        def narrow(ops):
            return Circuit(len(qa),
                           [Gate(g.name, g.params,
                                 tuple(remap[q] for q in g.qubits))
                            for g in ops])

        if not check_equivalent(narrow(sa_ops), narrow(sb_ops), tol):
            # inconclusive: this segmentation does not align the circuits —
            # never claim inequivalence from a misaligned cut
            return {"equivalent": None, "tier": 4,
                    "method": "compositional_segments",
                    "segments": len(segs_a), "undecided_segment": i,
                    "segment_qubits": qa, "global_phase_ignored": True}
    return {"equivalent": True, "tier": 4,
            "method": "compositional_segments",
            "segments": len(segs_a), "global_phase_ignored": True}
