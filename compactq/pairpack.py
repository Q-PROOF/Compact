"""Pair-locality packing — exact reordering that enlarges KAK windows.

Gates on disjoint wires commute exactly, so any topological order that
keeps each wire's sequence intact preserves the unitary.  `pair_pack`
chooses, at every step, the ready gate (all earlier gates on its wires
already emitted) that continues the current qubit pair when possible:
same-pair 2q gates pack into contiguous windows, so `kak.collect_blocks`
sees larger blocks and KAK resynthesis gets more to work with.

Exactness needs no proof beyond the disjoint-wire commutation rule and
per-wire order preservation — both are property-tested against the
dense prover.
"""
from __future__ import annotations

from .circuit import Circuit, Gate

__all__ = ["pair_pack"]


def _pair_of(g: Gate):
    """Sorted wire pair for a 2q gate, else None."""
    if len(g.qubits) == 2:
        return (min(g.qubits), max(g.qubits))
    return None


def pair_pack(circ: Circuit) -> Circuit:
    """Greedy pair-local topological reordering (exact).

    Gates are emitted in a stable greedy order: among the ready gates,
    prefer one whose qubit pair matches the previously emitted 2q gate;
    among ties, the earliest original position wins.  1q gates always
    yield to matching-pair 2q gates so blocks stay contiguous.
    """
    if len(circ.ops) < 3:
        return circ
    n_qubits = circ.num_qubits
    # per-wire queues (original order preserved per wire)
    queues: dict = {}
    for g in circ.ops:
        for q in g.qubits:
            queues.setdefault(q, []).append(g)

    emitted: set = set()
    next_idx: dict = {}          # gate object -> its per-wire position
    out: list = []
    last_pair = None
    total = len(circ.ops)

    # ready check: gate is ready when it is the head of every wire queue
    # it appears on
    def ready(g) -> bool:
        return all(queues[q][0] is g for q in g.qubits)

    while len(out) < total:
        best = None
        best_pair = None
        # scan wires in fixed order; heads are the candidate set
        for q in sorted(queues):
            if not queues[q]:
                continue
            g = queues[q][0]
            if not ready(g):
                continue
            p = _pair_of(g)
            if best is None:
                best, best_pair = g, p
                continue
            # prefer continuing the current 2q pair
            if last_pair is not None and p == last_pair:
                if best_pair != last_pair:
                    best, best_pair = g, p
                    continue
                # tie within the pair: keep the earlier wire's head
            elif best_pair == last_pair and p != last_pair:
                continue
            # among non-matching candidates prefer 2q gates (they unblock
            # pair sequences) then the lowest wire
            if best_pair != last_pair:
                if p is not None and best_pair is None:
                    best, best_pair = g, p
        if best is None:               # cannot happen on valid circuits
            break
        for q in best.qubits:
            queues[q].pop(0)
        out.append(best)
        emitted.add(id(best))
        p = _pair_of(best)
        if p is not None:
            last_pair = p

    return Circuit(n_qubits, out)
