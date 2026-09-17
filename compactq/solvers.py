"""Algorithmic solver layer: MaxCut-QAOA end to end.

The Fire-Opal-style product surface above the circuit level: given a
graph, build the QAOA ansatz, optimize the (gamma, beta) parameters
against the built-in exact simulator, execute the best circuit through
the suppression stack (or raw), and score the result against the
brute-force optimum.

Honest scope: exact-statevector parameter optimization is practical to
n <= 12; the classical scan is a coarse grid + local refinement (not a
high-performance optimizer).  What the layer demonstrates and gates is
the full path graph -> compiled QAOA -> suppression -> scored result.
"""
from __future__ import annotations

import itertools
import math

from .circuit import Circuit, Gate
from .simulate import exact_probabilities


def qaoa_circuit(edges, n, gamma, beta, p: int = 1) -> Circuit:
    """MaxCut-QAOA circuit on `n` wires (vertices 0..n-1): p rounds of
    [cost: CX RZ CX per edge] + [mixer: RX per wire].  Weights are 1
    and the (gamma, beta) pair repeats each round (linear schedule)."""
    ops = [Gate('h', (), (w,)) for w in range(n)]
    gs, bs = gamma, beta
    for _ in range(max(1, p)):
        for (a, b) in edges:
            ops += [Gate('cx', (), (a, b)),
                    Gate('rz', (2 * gs,), (b,)),
                    Gate('cx', (), (a, b))]
        for w in range(n):
            ops.append(Gate('rx', (2 * bs,), (w,)))
        gs, bs = gs * 1.0, bs * 1.0
    return Circuit(n, ops)


def cut_value(bitstr: str, edges) -> int:
    return sum(1 for (a, b) in edges if bitstr[a] != bitstr[b])


def brute_force_maxcut(edges, n):
    """(optimal_cut, optimal_bitstring) by exhaustive search (n <= 20)."""
    best, bs = -1, None
    for mask in range(1 << n):
        s = format(mask, "0%db" % n)[::-1]
        v = cut_value(s, edges)
        if v > best:
            best, bs = v, s
    return best, bs


def _score_probs(probs: dict, edges) -> float:
    tot = sum(probs.values()) or 1
    return sum(probs.get(s, 0) * cut_value(s, edges)
               for s in probs) / tot


def maxcut_qaoa(edges, *, p: int = 1, noise=None, shots: int = 4000,
                seed: int = 0, grid: int = 12, suppress: bool = True):
    """Solve MaxCut on `edges` with QAOA.

    Parameter loop: grid scan over (gamma, beta) on the built-in exact
    simulator (noiseless - the standard QAOA protocol), then execution
    of the best parameter set through the suppression stack when
    `suppress=True` (variant 0 circuit, pooled variants off; the
    measured distribution comes from the calibrated simulator).

    Returns best parameters, the scored distribution (mitigated,
    suppression applied when requested), the brute-force optimum and
    the approximation ratio."""
    if not edges:
        raise ValueError("graph has no edges")
    n = max(max(a, b) for (a, b) in edges) + 1
    best_cut, best_str = brute_force_maxcut(edges, n)

    best_params, best_score, best_probs = None, -1.0, None
    gammas = [math.pi * (i + 1) / (grid + 1) for i in range(grid)]
    betas = [math.pi * (i + 1) / (2 * (grid + 1)) for i in range(grid)]
    for g in gammas:
        for b in betas:
            probs = exact_probabilities(qaoa_circuit(edges, n, g, b, p))
            sc = _score_probs(probs, edges)
            if sc > best_score:
                best_params, best_score, best_probs = (g, b), sc, probs

    # execute the winning parameter set through suppression when asked
    if suppress:
        from .suppress import suppress_execute
        circ = qaoa_circuit(edges, n, *best_params, p)
        res = suppress_execute(circ, noise, seed=seed, variants=2,
                               shots=shots, mitigate=True)
        meas = res['probabilities']
    else:
        from .simulate import simulate_counts
        meas = simulate_counts(
            qaoa_circuit(edges, n, *best_params, p), noise, shots=shots,
            seed=seed)

    measured_cut = _score_probs(meas, edges)
    top = max(meas, key=meas.get)
    return {"edges": list(edges), "n": n, "p": p,
            "best_params": best_params,
            "ideal_expectation": best_score,
            "measured_expectation": measured_cut,
            "brute_force_optimum": best_cut,
            "best_bitstring": top,
            "approximation_ratio": measured_cut / best_cut,
            "optimal_found": cut_value(top, edges) == best_cut}
