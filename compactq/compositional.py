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


def _prefix_signatures(ops):
    """frozenset of wires touched by ops[:i] for every prefix cut i."""
    sigs = [frozenset()]
    cur = set()
    for g in ops:
        cur.update(g.qubits)
        sigs.append(frozenset(cur))
    return sigs


def verify_windows(a: Circuit, b: Circuit, max_width: int = 8,
                   max_checks: int = 2000, max_windows: int = 512,
                   tol: float = 1e-7):
    """T4.2 sliding-window compositional proof for CONNECTED circuits.

    Soundness: if the pair can be cut into matched windows (a's gates
    [lo..i) vs b's [lo'..j)) such that every window pair is equivalent
    up to global phase on the window's own wires, then the products are
    equivalent — U_a = S·W, U_b = S'·W', W ≡ W' reduces the question to
    the suffixes, recursively.  Every window is dense-verified on its
    OWN active wires (≤ `max_width`), so no 2^n unitary is ever built
    for the full circuit.  One-directional like verify_segmented:
    matched windows prove equivalence; an exhausted search DECLINES
    (None) — it is not evidence of inequivalence.

    The search enumerates candidate cut pairs whose windows have equal
    active wire sets, widest first, then proportional to the gate-count
    ratio; dead ends are memoized.  Bounded by `max_checks` dense
    comparisons and `max_windows`; returns None when the bounds or the
    width limit are exceeded, or when no alignment exists (callers fall
    through to the next prover).
    """
    from .equivalence import check_equivalent
    import sys as _sys

    if a.num_qubits != b.num_qubits:
        return None
    la, lb = len(a.ops), len(b.ops)
    if la == 0 and lb == 0:
        return {"equivalent": True, "tier": 4, "method": "sliding_windows",
                "windows": 0, "global_phase_ignored": True}

    # the DFS is deep (one frame per matched window); keep headroom and
    # restore the previous limit on the way out
    old_limit = _sys.getrecursionlimit()
    _sys.setrecursionlimit(max(old_limit, 10000))
    try:
        return _verify_windows_impl(a, b, la, lb, max_width, max_checks,
                                    max_windows, tol, check_equivalent)
    finally:
        _sys.setrecursionlimit(old_limit)


def _verify_windows_impl(a, b, la, lb, max_width, max_checks, max_windows,
                         tol, check_equivalent):
    checks = [0]
    failed = set()

    def window_ok(lo_a, i, lo_b, j):
        wa = sorted({q for g in a.ops[lo_a:i] for q in g.qubits})
        wb = sorted({q for g in b.ops[lo_b:j] for q in g.qubits})
        if not wa or wa != wb or len(wa) > max_width:
            return False
        remap = {q: x for x, q in enumerate(wa)}

        def narrow(ops, lo, hi):
            return Circuit(len(wa),
                           [Gate(g.name, g.params,
                                 tuple(remap[q] for q in g.qubits))
                            for g in ops[lo:hi]])

        checks[0] += 1
        return check_equivalent(narrow(a.ops, lo_a, i),
                                narrow(b.ops, lo_b, j), tol)

    def dfs(lo_a, lo_b):
        if (lo_a, lo_b) in failed:
            return None
        if lo_a == la and lo_b == lb:
            return 0
        # enumerate windows from lo_a: grow one gate at a time; the
        # window's ACTIVE set is the union of its own gates' wires
        # (gates may re-touch wires seen before earlier cuts — that is
        # fine; what is bounded is the window's own distinct-wire count)
        cands = []
        i = lo_a
        wa = set()
        while i < la:
            wa.update(a.ops[i].qubits)
            i += 1
            if len(wa) > max_width:
                break
            if lo_b == lb:
                continue
            j = lo_b
            wb = set()
            while j < lb:
                wb.update(b.ops[j].qubits)
                j += 1
                if len(wb) > max_width:
                    break
                if wb == wa:
                    cands.append((i + j, i, j))
        # prefer LARGE windows (shallower recursion, fewer dense checks),
        # then proportional cuts among comparable sizes: compilers
        # rewrite locally, so the matching cut in b sits near the global
        # gate-count ratio (both sides in window-relative positions)
        ratio = (lb - lo_b) / max(1, (la - lo_a))
        cands.sort(key=lambda t: (-(t[1] - lo_a + t[2] - lo_b),
                                  abs((t[2] - lo_b)
                                      - ratio * (t[1] - lo_a))))
        tried = 0
        for _score, i, j in cands:
            if tried >= 128 or checks[0] >= max_checks:
                break
            tried += 1
            if not window_ok(lo_a, i, lo_b, j):
                continue
            sub = dfs(i, j)
            if sub is not None:
                return sub + 1
        failed.add((lo_a, lo_b))
        if len(failed) > max_windows * 8:
            return None
        return None

    windows = dfs(0, 0)
    if windows is None:
        return None
    return {"equivalent": True, "tier": 4, "method": "sliding_windows",
            "windows": windows, "global_phase_ignored": True}
