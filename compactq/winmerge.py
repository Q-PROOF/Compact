"""Unitary-verified window merging.

For two same-pair KAK windows separated by gap gates, a gap gate g is moved
past the second window iff it *numerically* commutes with it (4x4 window
unitary comparison - no hand rules).  When every gap op moves, the windows
fuse and KAK can re-synthesize their product.  Whole-circuit fidelity
verified for <= 6 qubits.
"""
from __future__ import annotations

from .circuit import Circuit, Gate
from .equivalence import check_equivalent, unitary


def _window_commutes_with_op(w2_ops: list[Gate], g: Gate, remap: dict) -> bool:
    """True iff g commutes with the product of w2_ops (exact, 4x4).

    All gates are re-indexed through `remap` (original wires -> 0/1).
    """
    rem = lambda q: remap[q]
    w2r = [Gate(x.name, x.params, tuple(rem(q) for q in x.qubits)) for x in w2_ops]
    gr = Gate(g.name, g.params, tuple(rem(q) for q in g.qubits))
    wc = Circuit(2, w2r)
    u_w2 = unitary(wc)
    gc = Circuit(2, [gr])
    u_g = unitary(gc)
    # g after w2: u_g @ u_w2 ; w2 after g: u_w2 @ u_g
    left = [[sum(u_g[i][k] * u_w2[k][j] for k in range(4)) for j in range(4)]
            for i in range(4)]
    right = [[sum(u_w2[i][k] * u_g[k][j] for k in range(4)) for j in range(4)]
             for i in range(4)]
    tr = sum(left[i][j].conjugate() * right[i][j] for i in range(4) for j in range(4))
    return abs(tr) / 4 > 1 - 1e-9


def merge_by_unitary(circ: Circuit, max_steps: int = 32) -> Circuit:
    """Merge same-pair KAK windows across gap gates that numerically
    commute with the later window.

    One verified move at a time: find the first same-pair block pair whose
    gap ops ALL commute (4x4 unitary check) with the later window, splice
    them after that window, re-verify the whole circuit, and repeat.
    """
    if circ.num_qubits > 6:
        return circ
    cur = circ
    for _ in range(max_steps):
        ops = list(cur.ops)
        blocks = _collect_pair_blocks(ops)
        if len(blocks) < 2:
            break
        by_pair = {}
        for blk in blocks:
            by_pair.setdefault((blk[2], blk[3]), []).append(blk)
        applied = False
        for pair, blks in by_pair.items():
            for w1, w2 in zip(blks, blks[1:]):
                _, e1, _, _ = w1
                s2, e2, pa, pb = w2
                gap = list(range(e1, s2))
                if not gap or len(gap) > 64:
                    continue
                w2_ops = ops[s2:e2]
                remap = {pa: 0, pb: 1}
                movable = []
                ok = True
                for gpos in gap:
                    g = ops[gpos]
                    if set(g.qubits) != {pa, pb}:
                        ok = False
                        break
                    try:
                        if not _window_commutes_with_op(w2_ops, g, remap):
                            ok = False
                            break
                    except Exception:
                        ok = False
                        break
                    movable.append(gpos)
                if not (ok and movable):
                    continue
                # splice: [.. W1 ..][gap -> moved after W2][W2][rest]
                new_ops = ops[:e1] + ops[s2:e2]
                for gpos in gap:
                    if gpos not in movable:
                        new_ops.append(ops[gpos])
                new_ops.extend(ops[e2:])
                new_c = Circuit(circ.num_qubits, new_ops)
                if check_equivalent(cur, new_c):
                    cur = new_c
                    applied = True
                break
            if applied:
                break
        if not applied:
            break
    return cur


def _collect_pair_blocks(ops):
    """(start, end, a, b) maximal same-pair windows (with 1q on pair)."""
    n_ops = len(ops)
    used = [False] * n_ops
    blocks = []
    i = 0
    while i < n_ops:
        if not used[i] and len(ops[i].qubits) == 2:
            a, b = ops[i].qubits
            wires = {a, b}
            j = i
            while j < n_ops and not used[j] and set(ops[j].qubits) <= wires:
                j += 1
            start = i
            while start > 0 and not used[start - 1] \
                    and len(ops[start - 1].qubits) == 1 \
                    and ops[start - 1].qubits[0] in wires:
                start -= 1
            for k in range(start, j):
                used[k] = True
            blocks.append((start, j, min(a, b), max(a, b)))
            i = j
        else:
            i += 1
    return blocks
