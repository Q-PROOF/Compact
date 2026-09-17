"""Approximate optimization: trade a bounded, quantified amount of fidelity
for fewer two-qubit gates.

Two-qubit blocks are re-synthesized into the cheapest CX class whose average
gate fidelity |Tr(U_ref^dag U)|/d (d = 4 for a two-qubit block) is at least
`min_fidelity`.  Each synthesized block carries that per-block guarantee from
the exact KAK fidelity table; the whole-circuit infidelity is therefore
bounded by (number of approximated blocks) * (1 - min_fidelity).
"""
from __future__ import annotations

from .circuit import Circuit
from .search import optimize_search


def approximate(circ: Circuit, min_fidelity: float = 0.999) -> Circuit:
    """Approximately optimize `circ`, prioritizing 2-qubit gate count.

    Per-block average gate fidelity is guaranteed >= min_fidelity; total
    circuit infidelity is bounded by (#approximated blocks) * (1 -
    min_fidelity).  An exact baseline is always computed first and an
    approximate candidate is returned only when it is strictly smaller.
    """
    if not 0.0 < min_fidelity <= 1.0:
        raise ValueError("min_fidelity must be in (0, 1]")
    return optimize_search(circ, verify=None,
                           fidelity_tolerance=min_fidelity)
