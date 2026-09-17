"""Clifford Data Regression (CDR) - learnable error mitigation.

The near-Clifford trick: build K variants of the circuit in which each
non-Clifford rotation is SNAPPED to the nearest Clifford angle with
some probability.  Every variant is executed on the noisy target (or
the calibrated simulator) AND evaluated exactly on the noiseless
simulator.  The (noisy, exact) pairs span the noise response of the
observable near this circuit, so a simple linear fit

    E_exact ~ alpha * E_noisy + beta

learned on the variants corrects the noisy expectation of the ORIGINAL
circuit.  CDR is the strongest mitigation when the circuit is close to
Clifford - exactly the regime where our stabilizer machinery gives
exact training labels for free, with no device shots burned on the
exact side.

Zero-dependency.  Composes with the suppression stack (run the
variants through the same run_fn) and with ZNE (different correction
axis: CDR learns the noise response from structured variants, ZNE
scales the noise).
"""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate

# 1q gates that are already Clifford (no snapping needed)
_CLIFFORD_1Q = {"h", "x", "y", "z", "s", "sdg", "sx", "sxdg"}


def _is_clifford_angle(theta: float) -> bool:
    k = round(theta / (math.pi / 2))
    return abs(theta - k * math.pi / 2) < 1e-9


def _snap(theta: float) -> float:
    return round(theta / (math.pi / 2)) * (math.pi / 2)


def near_clifford_variants(circ: Circuit, k: int, rng: random.Random,
                           p: float = 0.65):
    """Build `k` near-Clifford variants: every non-Clifford rotation is
    snapped to its nearest Clifford angle independently with
    probability `p` (per gate, per variant), so the family spans the
    line between the Clifford ideal and the original circuit."""
    variants = []
    for _ in range(k):
        ops = []
        for g in circ.ops:
            if g.name in ("rz", "rx", "ry", "p") and g.params \
                    and not _is_clifford_angle(float(g.params[0])) \
                    and rng.random() < p:
                ops.append(Gate(g.name, (_snap(float(g.params[0])),),
                                g.qubits))
            elif g.name == "t" and rng.random() < p:
                ops.append(Gate("s", (), g.qubits))       # t ~ rz(pi/4) -> s
            elif g.name == "tdg" and rng.random() < p:
                ops.append(Gate("sdg", (), g.qubits))
            else:
                ops.append(g)
        variants.append(Circuit(circ.num_qubits, ops))
    return variants


def parity_expectation_from_counts(counts: dict, obs_wires, n: int) -> float:
    """E = sum_s (-1)^{popcount(s restricted to obs_wires)} p_s."""
    mask = [0] * n
    for w in obs_wires:
        mask[w] = 1
    tot = sum(counts.values()) or 1
    acc = 0.0
    for s, c in counts.items():
        par = sum(int(s[w]) * mask[w] for w in range(n)) & 1
        acc += (1 if par == 0 else -1) * c
    return acc / tot


def parity_expectation_exact(probs: dict, obs_wires, n: int) -> float:
    mask = [0] * n
    for w in obs_wires:
        mask[w] = 1
    acc = 0.0
    for s, p in probs.items():
        par = sum(int(s[w]) * mask[w] for w in range(n)) & 1
        acc += (1 if par == 0 else -1) * p
    return acc


def _linear_fit(xs, ys):
    """Least-squares y ~ alpha x + beta (2-parameter normal equations)."""
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx < 1e-15:
        return 0.0, my  # degenerate: predict the mean
    alpha = sxy / sxx
    return alpha, my - alpha * mx


def cdr_execute(circ: Circuit, obs_wires, noise=None, run_fn=None, *,
                training_size: int = 6, shots: int = 4000, seed: int = 0,
                snap_probability: float = 0.65,
                eps_coh: float = 0.0) -> dict:
    """One-call CDR for the Z-parity observable on `obs_wires`.

    Training: `training_size` near-Clifford variants are executed
    NOISILY (through `run_fn` on hardware, or the calibrated built-in
    simulator) and evaluated EXACTLY on the noiseless simulator - the
    exact side never burns device shots.  A linear (noisy -> exact)
    map is fit and applied to the original circuit's noisy expectation.

    Returns the corrected expectation, the unmitigated value, the fit
    parameters and the per-variant training table."""
    from .simulate import simulate_counts, exact_probabilities
    noise = noise  # may be None only when run_fn is provided

    rng = random.Random(seed)
    variants = near_clifford_variants(circ, training_size, rng,
                                      p=snap_probability)
    n = circ.num_qubits

    def noisy_counts(c):
        if run_fn is not None:
            return run_fn(c, seed=rng.randrange(1 << 30), shots=shots)
        return simulate_counts(c, noise, shots=shots, seed=rng.randrange(1 << 30),
                               eps_coh=eps_coh)

    xs, ys = [], []
    for v in variants:
        e_noisy = parity_expectation_from_counts(noisy_counts(v), obs_wires, n)
        e_exact = parity_expectation_exact(exact_probabilities(v),
                                           obs_wires, n)
        xs.append(e_noisy)
        ys.append(e_exact)

    alpha, beta = _linear_fit(xs, ys)

    e0_noisy = parity_expectation_from_counts(noisy_counts(circ), obs_wires, n)
    e0_exact = parity_expectation_exact(exact_probabilities(circ),
                                        obs_wires, n)
    corrected = alpha * e0_noisy + beta
    corrected = max(-1.0, min(1.0, corrected))  # physical range clamp
    return {"expectation": corrected,
            "unmitigated": e0_noisy,
            "ideal_reference": e0_exact,
            "alpha": alpha,
            "beta": beta,
            "training": [{"noisy": x, "exact": y} for x, y in zip(xs, ys)],
            "snap_probability": snap_probability,
            "training_size": training_size}
