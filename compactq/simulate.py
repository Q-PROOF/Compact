"""Zero-dependency stochastic-wavefunction noise simulator.

One trajectory applies the circuit to a pure statevector and samples the
jump processes; ``simulate_counts`` repeats for `shots` and bins the
measured bitstrings.  Error classes mirror the density-matrix study in
``scripts/suppress_bench.py`` so results transfer:

  - stochastic gate depolarizing (per wire, from NoiseModel infidelity)
  - T1 relaxation and dephasing jumps during idle windows (ASAP schedule
    from gate durations, so DD pulses genuinely consume idle time)
  - quasi-static Z drift: one rate per wire per call, sign-flipped by
    every X/Y pulse - the exact mechanism dynamical decoupling exploits
  - coherent overrotation per 2q gate: axis drawn once per call
    (quasi-static session - this is what randomized compiling averages
    over VARIANTS, so raw and twirled submissions differ systematically)
  - readout confusion on the measured bitstring

Pure Python (no numpy in the core): practical to n <= 14, shots in the
thousands.  Bitstrings use compactq's little-endian convention: char i
of the string is wire i.
"""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate
from .linalg import gate_matrix
from .noise import NoiseModel, default_model
from .equivalence import unitary

_I2 = (1 + 0j, 0j, 0j, 1 + 0j)


def _mat2(g: Gate):
    m = gate_matrix(g)
    return tuple(complex(x) for x in m)


def _mat4(g: Gate):
    rows = unitary(Circuit(2, [Gate(g.name, g.params, (0, 1))]))
    return tuple(complex(x) for row in rows for x in row)


def _apply_1q(state, n, w, m):
    step = 1 << w
    a, b, c, d = m
    for base in range(0, len(state), step * 2):
        for j in range(base, base + step):
            x, y = state[j], state[j + step]
            state[j] = a * x + b * y
            state[j + step] = c * x + d * y


def _apply_2q(state, n, w0, w1, m):
    s0, s1 = 1 << w0, 1 << w1
    idx = [0] * 4
    out = [0j] * 4
    for j in range(len(state)):
        b0 = (j >> w0) & 1
        b1 = (j >> w1) & 1
        if b0 or b1:
            continue  # process only the sub-space base
        idx[0] = j
        idx[1] = j | s0
        idx[2] = j | s1
        idx[3] = j | s0 | s1
        v = (state[idx[0]], state[idx[1]], state[idx[2]], state[idx[3]])
        for r in range(4):
            acc = 0j
            row = m[r * 4:(r + 1) * 4]
            for cc in range(4):
                acc += row[cc] * v[cc]
            out[r] = acc
        for r in range(4):
            state[idx[r]] = out[r]


def _apply_pauli(state, n, w, name):
    if name == "x":
        _apply_1q(state, n, w, (0 + 0j, 1 + 0j, 1 + 0j, 0 + 0j))
    elif name == "y":
        _apply_1q(state, n, w, (0j, -1j, 1j, 0j))
    elif name == "z":
        step = 1 << w
        for base in range(0, len(state), step * 2):
            for j in range(base, base + step):
                state[j + step] = -state[j + step]


def _rz_phase(state, n, w, theta):
    half = theta / 2.0
    epos = complex(math.cos(half), -math.sin(half))
    eneg = complex(math.cos(half), math.sin(half))
    step = 1 << w
    for base in range(0, len(state), step * 2):
        for j in range(base, base + step):
            state[j] *= epos
            state[j + step] *= eneg


def _prob_of(state):
    return [abs(x) * abs(x) for x in state]


def sample_outcome(state, rng):
    p = _prob_of(state)
    total = sum(p)
    if total <= 0:
        return 0
    r = rng.random() * total
    acc = 0.0
    for j, pj in enumerate(p):
        acc += pj
        if r < acc:
            return j
    return len(p) - 1


def run_shot(circ: Circuit, noise: NoiseModel, rng: random.Random,
             eps_coh: float = 0.0, axes=None) -> int:
    """One noisy trajectory; returns the raw (pre-readout) outcome index."""
    n = circ.num_qubits
    state = [0j] * (1 << n)
    state[0] = 1.0 + 0j
    t = [0.0] * n
    dsign = [1] * n
    drift = []
    for w in range(n):
        s = noise.drift(w)
        drift.append(rng.gauss(0.0, s) if s > 0 else 0.0)
    if eps_coh > 0:
        if axes is None:
            axes = [rng.choice(("x", "y", "z")) for _ in range(256)]
    i2q = 0
    cache2 = {}

    def idle(w, dt):
        if dt <= 0:
            return
        if drift[w] != 0.0:
            _rz_phase(state, n, w, drift[w] * dt * dsign[w])
        t1_ns = noise.t1(w) * 1000.0
        gamma = 1.0 - math.exp(-dt / t1_ns)
        if gamma > 0 and rng.random() < gamma:
            # relaxation jump: amplitude-damping K1 maps |1> -> |0>
            step = 1 << w
            for j in range(len(state)):
                if (j >> w) & 1:
                    state[j & ~(1 << w)] += state[j]
                    state[j] = 0j
            norm = 0.0
            for x in state:
                norm += abs(x) * abs(x)
            if norm > 0:
                nrm = math.sqrt(norm)
                for j in range(len(state)):
                    state[j] /= nrm
        t2_ns = noise.t2(w) * 1000.0
        tphi_inv = max(0.0, 1.0 / t2_ns - 0.5 / t1_ns)
        lam = 1.0 - math.exp(-dt * tphi_inv)
        if lam > 0 and rng.random() < lam:
            _apply_pauli(state, n, w, "z")

    for g in circ.ops:
        dur = noise.duration_of(g)
        start_t = max(t[w] for w in g.qubits) if g.qubits else 0.0
        for w in g.qubits:
            idle(w, start_t - t[w])
        if len(g.qubits) == 1:
            _apply_1q(state, n, g.qubits[0], _mat2(g))
        else:
            m = cache2.get((g.name, g.params))
            if m is None:
                m = _mat4(g)
                cache2[(g.name, g.params)] = m
            _apply_2q(state, n, g.qubits[0], g.qubits[1], m)
        # stochastic depolarizing per involved wire
        p = noise.gate_error(g)
        for w in g.qubits:
            if rng.random() < p / max(1, len(g.qubits)):
                _apply_pauli(state, n, w,
                             rng.choice(("x", "y", "z")))
        if len(g.qubits) == 2 and eps_coh > 0:
            ax = axes[i2q % len(axes)]
            i2q += 1
            c, s = math.cos(eps_coh / 2), math.sin(eps_coh / 2)
            if ax == "z":
                _rz_phase(state, n, g.qubits[0], eps_coh)
            elif ax == "x":   # Rx(theta)
                _apply_1q(state, n, g.qubits[0],
                          (c + 0j, -1j * s, -1j * s, c + 0j))
            else:             # Ry(theta)
                _apply_1q(state, n, g.qubits[0],
                          (c + 0j, -s + 0j, s + 0j, c + 0j))
        for w in g.qubits:
            t[w] = start_t + dur
            if g.name in ("x", "y"):
                dsign[w] = -dsign[w]
    total = max(t) if t else 0.0
    for w in range(n):
        idle(w, total - t[w])
    return sample_outcome(state, rng)


def simulate_counts(circ: Circuit, noise: NoiseModel | None = None,
                    shots: int = 1024, seed: int = 0,
                    eps_coh: float = 0.0) -> dict:
    """Sample `shots` noisy trajectories of `circ`; returns
    bitstring -> count (char i = wire i, little-endian)."""
    noise = noise or default_model(circ.num_qubits)
    rng = random.Random(seed)
    axes = ([rng.choice(("x", "y", "z")) for _ in range(256)]
            if eps_coh > 0 else None)
    counts: dict = {}
    for _ in range(shots):
        j = run_shot(circ, noise, rng, eps_coh=eps_coh, axes=axes)
        bits = "".join((j >> w) & 1 and "1" or "0"
                       for w in range(circ.num_qubits))
        counts[bits] = counts.get(bits, 0) + 1
    return counts


def statevector(circ: Circuit):
    """Exact noiseless statevector of `circ` (list of 2^n complex
    amplitudes, index j with wire w = bit w).  Pure Python; practical
    to n <= 14.  Deterministic - no sampling involved."""
    n = circ.num_qubits
    state = [0j] * (1 << n)
    state[0] = 1.0 + 0j
    cache2 = {}
    for g in circ.ops:
        if len(g.qubits) == 1:
            _apply_1q(state, n, g.qubits[0], _mat2(g))
        else:
            key = (g.name, g.params)
            m = cache2.get(key)
            if m is None:
                m = _mat4(g)
                cache2[key] = m
            _apply_2q(state, n, g.qubits[0], g.qubits[1], m)
    return state


def exact_probabilities(circ: Circuit) -> dict:
    """Exact ideal outcome distribution of `circ` (no sampling):
    bitstring -> probability (char i = wire i), zero-amplitude entries
    dropped.  Practical to n <= 14."""
    n = circ.num_qubits
    sv = statevector(circ)
    out = {}
    for j, amp in enumerate(sv):
        p = (amp.real * amp.real + amp.imag * amp.imag)
        if p > 1e-15:
            out[format(j, "0%db" % n)[::-1]] = p
    return out
