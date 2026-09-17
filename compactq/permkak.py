"""Permutation-absorbing KAK block resynthesis.

Each maximal same-pair KAK block is synthesized with a free content swap:
the block may act on its wires in either orientation, whichever synthesizes
to fewer 2q gates (the swapped orientation can land the block's Weyl
coordinates in a cheaper class).  Choosing per-block minima lets qubit
content hop between wires across the circuit - exactly the freedom Qiskit's
L3 exploits via permutation elision - and the residual wire permutation is
paid once, at the end, as explicit SWAP gates.

Exactness: the emitted circuit is re-verified against the input with the
dense prover before it can replace anything (compactq's proof net).
"""
from __future__ import annotations

from .circuit import Circuit, Gate


def _relabel(gates, mapping, n):
    return [Gate(g.name, g.params, tuple(mapping.get(w, w) for w in g.qubits))
            for g in gates]


def _permutation_swaps(pos):
    """SWAP transpositions (physical wires) restoring content to home wires.

    `pos[w]` = physical wire currently holding original wire w's content.
    """
    pos = list(pos)
    swaps = []
    for w in range(len(pos)):
        while pos[w] != w:
            cur = pos[w]
            pos[w], pos[cur] = pos[cur], pos[w]
            swaps.append((w, cur))
    return swaps


def permutation_kak_pass(circ: Circuit, verify: bool = True):
    """Returns an exact Circuit with permutation-absorbed KAK blocks, or None
    when nothing improved (or the prover rejected the rebuild)."""
    from .kak import collect_blocks
    blocks = collect_blocks(circ)
    if not blocks:
        return None
    n = circ.num_qubits
    ops = list(circ.ops)
    pos = list(range(n))          # original wire w -> physical wire
    out = []
    prev = 0
    changed = False

    for (s, e, a, b) in blocks:
        out.extend(_relabel(ops[prev:s], {w: pos[w] for w in range(n)}, n))
        prev = e
        pa, pb = pos[a], pos[b]
        id_map = {w: pos[w] for w in range(n)}
        x_map = dict(id_map)
        x_map[a], x_map[b] = pb, pa

        seg_id = _relabel(ops[s:e], id_map, n)
        seg_x = _relabel(ops[s:e], x_map, n)
        bc_id = Circuit(n, seg_id)
        bc_x = Circuit(n, seg_x)
        cx_id = bc_id.two_qubit_count()
        cx_x = bc_x.two_qubit_count()
        if cx_x < cx_id:
            out.extend(seg_x)
            changed = True
            # content of a and b crossed: update the ledger
            for w in range(n):
                if pos[w] == pa:
                    pos[w] = pb
                elif pos[w] == pb:
                    pos[w] = pa
        else:
            out.extend(seg_id)
    out.extend(_relabel(ops[prev:], {w: pos[w] for w in range(n)}, n))

    # residual permutation: pay it once, as SWAPs
    swaps = _permutation_swaps(pos)
    if swaps:
        changed = True
    for a, b in swaps:
        out.extend([Gate("cx", (), (a, b)), Gate("cx", (), (b, a)),
                    Gate("cx", (), (a, b))])

    if not changed:
        return None
    rebuilt = Circuit(n, out)
    if rebuilt.two_qubit_count() >= circ.two_qubit_count():
        return None
    if verify:
        from .equivalence import check_equivalent
        if not check_equivalent(circ, rebuilt):
            return None
    return rebuilt
