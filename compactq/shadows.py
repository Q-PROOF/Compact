"""Classical shadows: many observables from one measurement budget.

Protocol (Huang-Kuang-Preskill style, single-qubit Pauli ensemble):
each snapshot applies a uniformly random single-qubit basis change per
wire, measures in Z, and stores the outcome.  The estimator for a
Z-parity observable on wire set W: a snapshot contributes only when
every wire of W was measured in the Z basis ("matched basis" rule),
scaled by the randomized-measurement factor 3^|W|; the reported value
is a median-of-means over snapshot groups for heavy-tail robustness.

Zero dependencies.  Works with any run_fn (hardware adapters) or the
built-in simulator.
"""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate

# basis-change gates per measurement basis
_BASIS_GATE = {"z": (), "x": ("h",), "y": (("sdg"),)}


def _basis_rotation(basis: str) -> tuple:
    """Gates mapping the measurement basis onto the Z basis."""
    if basis == "x":
        return ("h",)
    if basis == "y":
        return ("sdg", "h")   # measure Y: rotate Y -> Z
    return ()


def shadow_snapshots(circ: Circuit, n_snapshots: int, shots_per: int,
                     noise=None, rng: random.Random | None = None,
                     seed: int = 0):
    """Run `n_snapshots` randomized-basis measurement rounds.  Returns
    snapshots as (basis_tuple, {bitstring: count}) pairs.  `basis_tuple`
    is the per-wire measurement basis ('x'/'y'/'z')."""
    from .simulate import simulate_counts
    rng = rng or random.Random(seed)
    n = circ.num_qubits
    snaps = []
    for _ in range(n_snapshots):
        bases = tuple(rng.choice(("x", "y", "z")) for _ in range(n))
        ops: list = []
        for w in range(n):
            for gname in _basis_rotation(bases[w]):
                ops.append(Gate(gname, (), (w,)))
        probe = Circuit(n, ops + list(circ.ops))
        counts = simulate_counts(probe, noise, shots=shots_per,
                                 seed=rng.randrange(1 << 30))
        snaps.append((bases, counts))
    return snaps


def shadow_estimate_parity(snaps, obs_wires, n: int,
                           groups: int = 5) -> dict:
    """Estimate the Z-parity observable on `obs_wires` from snapshots.

    Snapshot contribution (matched wires only, others contribute 0):
    product over observed wires of (3 * outcome_sign); the UNCONDITIONAL
    mean over all snapshots is the unbiased classical-shadow estimate.
    Median-of-means over `groups` gives a heavy-tail-robust variant."""
    k = len(obs_wires)
    contribs = []
    for (bases, counts) in snaps:
        tot = sum(counts.values()) or 1
        if any(bases[w] != 'z' for w in obs_wires):
            contribs.append(0.0)
            continue
        acc = 0.0
        for s, c in counts.items():
            par = sum(int(s[w]) for w in obs_wires) & 1
            acc += (1 if par == 0 else -1) * c / tot
        contribs.append((3 ** k) * acc)
    est_uncond = sum(contribs) / max(1, len(contribs))
    per = max(1, len(contribs) // groups)
    means = []
    for g in range(groups):
        chunk = contribs[g * per:(g + 1) * per] if g < groups - 1 \
            else contribs[g * per:]
        if chunk:
            means.append(sum(chunk) / len(chunk))
    means.sort()
    median = means[len(means) // 2] if len(means) % 2 \
        else 0.5 * (means[len(means) // 2 - 1] + means[len(means) // 2])
    matched = sum(1 for c in contribs if c != 0.0)
    return {"expectation": est_uncond,
            "median_of_means": median,
            "matched_snapshots": matched,
            "total_snapshots": len(snaps)}
