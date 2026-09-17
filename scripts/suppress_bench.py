"""Error-suppression benchmark: Compact pipeline vs default execution.

A controlled noise simulator (numpy density matrix, n <= 5) implementing
the error classes commercial suppression pipelines target:

  - coherent overrotation (epsilon per 2q gate, applied immediately
    after the gate, i.e. between the gate and any following twirl
    correction - the physical reason randomized compiling works)
  - stochastic depolarizing noise per gate (from the noise model)
  - T1/T2 decoherence during idle windows (ASAP schedule, durations)
  - per-wire readout confusion on the final distribution

Pipelines compared (cumulative layers, mirroring the Mundada et al. 2023
protocol of "default execution vs suppressed execution"):

  raw          : the input circuit as-is
  optimized    : compactq.optimize_search (the compiler layer)
  +dd          : dynamical decoupling, benefit-gated per idle window
  +twirl       : Pauli twirling, K variants averaged (randomized compiling)
  +mem         : measurement error mitigation on the readout

Success metric per family: probability of the ideal outcome string after
readout.  The gate for the whole study: the full pipeline must beat raw
execution (factor > 1) in EVERY scenario.  This is a simulator-based,
fully reproducible measurement of deterministic suppression; hardware
runs need backend credentials and use the same harness functions.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from compactq import Circuit, Gate, optimize_search, benchmarks  # noqa: E402
from compactq.noise import NoiseModel  # noqa: E402
from compactq.suppress import (pauli_twirl, insert_dd, idle_windows,  # noqa: E402
                               expand_for_suppression)
from compactq.mitigate import mitigate_counts  # noqa: E402
from compactq.equivalence import unitary as unitary_of  # noqa: E402

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def rx(theta):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)


def full_gate(n, g):
    """Full 2^n unitary of a 1q/2q gate under compactq's little-endian
    wire convention (wire w = bit w)."""
    from compactq.linalg import gate_matrix
    dim = 1 << n
    U = np.zeros((dim, dim), dtype=complex)
    wires = list(g.qubits)
    if len(wires) == 1:
        m = np.asarray(gate_matrix(g), dtype=complex).reshape(2, 2)
        for j in range(dim):
            bit = (j >> wires[0]) & 1
            for r in range(2):
                k = (j & ~(1 << wires[0])) | (r << wires[0])
                U[k, j] += m[r, bit]
        return U
    a, b = wires
    # G4's LSB is its wire 0 = circuit wire a (little-endian throughout)
    G4 = np.asarray(unitary_of(Circuit(2, [Gate(g.name, g.params, (0, 1))])), dtype=complex)
    for j in range(dim):
        jb = ((j >> a) & 1) | (((j >> b) & 1) << 1)
        for r in range(4):
            k = (j & ~(1 << a) & ~(1 << b)) | ((r & 1) << a) | (((r >> 1) & 1) << b)
            U[k, j] += G4[r, jb]
    return U


def depol_1q(rho, w, n, p):
    if p <= 0:
        return rho
    Ux = full_gate(n, Gate("x", (), (w,)))
    Uy = full_gate(n, Gate("y", (), (w,)))
    Uz = full_gate(n, Gate("z", (), (w,)))
    return (1 - p) * rho + (p / 3.0) * (Ux @ rho @ Ux.conj().T
                                        + Uy @ rho @ Uy.conj().T
                                        + Uz @ rho @ Uz.conj().T)


def amp_damp_1q(rho, w, n, gamma):
    """Kraus amplitude damping on wire w."""
    if gamma <= 0:
        return rho
    dim = 1 << n
    K0 = np.zeros((dim, dim), dtype=complex)
    K1 = np.zeros((dim, dim), dtype=complex)
    for j in range(dim):
        if (j >> w) & 1 == 0:
            K0[j, j] = 1.0
            K1[j | (1 << w), j] = np.sqrt(gamma)
        else:
            K0[j, j] = np.sqrt(1 - gamma)
    return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T


def phase_damp_1q(rho, w, n, lam):
    if lam <= 0:
        return rho
    Uz = full_gate(n, Gate("z", (), (w,)))
    return (1 - lam) * rho + lam * (Uz @ rho @ Uz.conj().T)


def run_noisy(circ, noise, eps_coh, idle_decoherence=True,
              drift_sigma=None, seed=0):
    """Density matrix after noisy execution (pre-readout).

    Time-ordered: idle T1/T2 decay is applied per wire when its next gate
    arrives (and flushed at the end), so DD gates genuinely consume idle
    time.  Coherent overrotation fires immediately after each 2q gate -
    between the gate and any following twirl correction, which is exactly
    the mechanism randomized compiling exploits.  Quasi-static Z drift is
    drawn per wire from the noise model's drift_rate (or the drift_sigma
    override); each X/Y pulse flips the drift axis sign, which is exactly
    what XY4 decoupling exploits."""
    n = circ.num_qubits
    rho = np.zeros((1 << n, 1 << n), dtype=complex)
    rho[0, 0] = 1.0
    t = [0.0] * n

    rng = np.random.default_rng(seed)
    drift = {}
    for w in range(n):
        s = drift_sigma if drift_sigma is not None else noise.drift(w)
        drift[w] = rng.normal(0.0, s) if s > 0 else 0.0
    # coherent-error AXIS per 2q gate site, quasi-static within an
    # instance (calibration error is directional and varies per site).
    # Keyed by 2q-gate ordinal so every pipeline row sees the SAME axes:
    # rewrites (twirl/DD) add only 1q gates, never 2q ones.
    axes = [("rx", "ry", "rz")[int(rng.integers(0, 3))] for _ in range(256)]
    # sign of the drift axis on each wire: an X or Y pulse conjugates
    # Z -> -Z, so DD sequences flip the sign between segments and the
    # drift refocuses - this is exactly why XY4 works
    dsign = {w: 1 for w in range(n)}

    def decay(w, dt_ns):
        if dt_ns <= 0:
            return
        if idle_decoherence:
            t1_ns = noise.t1(w) * 1000.0
            t2_ns = noise.t2(w) * 1000.0
            gamma = 1 - np.exp(-dt_ns / t1_ns)
            tphi_inv = max(0.0, 1.0 / t2_ns - 0.5 / t1_ns)
            lam = (1 - np.exp(-dt_ns * tphi_inv)) / 2.0
            if gamma > 0:
                rho[:] = amp_damp_1q(rho, w, n, min(gamma, 0.999))
            if lam > 0:
                rho[:] = phase_damp_1q(rho, w, n, min(lam, 0.499))
        if drift[w] != 0.0:
            ang = drift[w] * dt_ns * dsign[w]
            Uz = full_gate(n, Gate("rz", (ang,), (w,)))
            rho[:] = Uz @ rho @ Uz.conj().T

    i2q = 0
    for g in circ.ops:
        dur = noise.duration_of(g)
        start_t = max(t[w] for w in g.qubits) if g.qubits else 0.0
        for w in g.qubits:
            decay(w, start_t - t[w])
        U = full_gate(n, g)
        rho[:] = U @ rho @ U.conj().T
        p = noise.gate_error(g)
        for w in g.qubits:
            rho[:] = depol_1q(rho, w, n, p / max(1, len(g.qubits)))
        if len(g.qubits) == 2 and eps_coh > 0:
            Ue = full_gate(n, Gate(axes[i2q], (eps_coh,), (g.qubits[0],)))
            rho[:] = Ue @ rho @ Ue.conj().T
            i2q += 1
        for w in g.qubits:
            t[w] = start_t + dur
            if g.name in ("x", "y"):
                dsign[w] = -dsign[w]  # pulse flips the drift axis
    total = max(t)
    for w in range(n):
        decay(w, total - t[w])
    return rho


def readout_apply(probs, noise, n):
    out = np.zeros_like(probs)
    for j in range(len(probs)):
        p = probs[j]
        if p == 0:
            continue
        for k in range(len(probs)):
            q = p
            ok = True
            for w in range(n):
                jb, kb = (j >> w) & 1, (k >> w) & 1
                p10, p01 = noise.readout_error(w)
                q *= (p10 if (jb, kb) == (0, 1) else
                      p01 if (jb, kb) == (1, 0) else
                      1 - (p10 if jb == 0 else p01))
            out[k] += q
    return out


def success(rho, good_mask, noise, mem=False):
    n = int(np.log2(rho.shape[0]))
    probs = np.real(np.diag(rho))
    probs = np.clip(probs, 0, None)
    probs /= probs.sum()
    rp = readout_apply(probs, noise, n)
    if mem:
        counts = {format(k, "0%db" % n)[::-1]: float(rp[k])
                  for k in range(len(rp)) if rp[k] > 0}
        # compactq mitigate expects char i = wire i; format gives MSB-first,
        # reverse for little-endian
        mit = mitigate_counts(counts, {w: noise.readout_error(w)
                                       for w in range(n)}, n, clip=True)
        vals = np.zeros(len(rp))
        for s, p in mit.items():
            vals[int(s[::-1], 2)] = p
        rp = vals
    return float(rp @ good_mask)


def good_mask_for(family, n, circ=None):
    """Success mask.  With `circ`, derived from the circuit's own ideal
    (zero-noise) distribution - robust to circuit-convention differences."""
    d = 1 << n
    m = np.zeros(d)
    if circ is not None:
        rho = np.zeros((d, d), dtype=complex)
        rho[0, 0] = 1.0
        for g in circ.ops:
            U = full_gate(n, g)
            rho = U @ rho @ U.conj().T
        diag = np.real(np.diag(rho))
        top = diag.max()
        m[diag > top - 1e-9] = 1.0
        return m
    if family == "ghz":
        m[0] = m[d - 1] = 1.0
    elif family == "bv":
        m[d - 1] = 1.0
    elif family == "qft":
        m[0] = 1.0
    return m


def build(family, n):
    if family == "ghz":
        ops = [Gate("h", (), (0,))]
        for j in range(n - 1):
            ops.append(Gate("cx", (), (j, j + 1)))
        return Circuit(n, ops)
    if family == "bv":
        # Z-oracle BV: ancilla in |1> (phase kickback needs no |-> here)
        ops = [Gate("x", (), (n - 1,))]
        ops += [Gate("h", (), (j,)) for j in range(n - 1)]
        for j in range(n - 1):
            ops.append(Gate("cz", (), (j, n - 1)))
        ops += [Gate("h", (), (j,)) for j in range(n - 1)]
        return Circuit(n, ops)
    if family == "qft":
        # X on wire 0 breaks the uniform |000> input symmetry: QFT|100..0>
        # is a phase-gradient state with a single deterministic outcome
        ops = [Gate("x", (), (0,))]
        return Circuit(n, ops + list(benchmarks.qft(n).ops))
    raise ValueError(family)




def scenario_model(n, scenario):
    """Per-mechanism stress models: each layer is measured against the
    error class it exists to suppress.  Drift (rms rad/ns) lives in the
    model itself so the benefit-gated DD pass can read the same physics
    the simulator applies."""
    m = NoiseModel(n)
    for w in range(n):
        m.t1_us[w] = 150.0
        m.t2_us[w] = 100.0
        m.readout[w] = (0.002, 0.002)
    m.gate_infidelity = {"1q": 0.0002, "cx": 0.001, "cz": 0.001,
                         "swap": 0.003, "cp": 0.001}
    drift = 0.0
    if scenario == "coherent":       # twirling's target
        m.gate_infidelity = {"1q": 0.0002, "cx": 0.0005, "cz": 0.0005,
                             "swap": 0.0015, "cp": 0.0005}
    elif scenario == "decoherence":  # DD's target: fast T1/T2, slow gates
        for w in range(n):
            m.t1_us[w] = 12.0
            m.t2_us[w] = 8.0
        m.durations_ns.update({"cx": 600.0, "cz": 600.0, "1q": 120.0})
        drift = 0.0012
    elif scenario == "readout":      # MEM's target
        for w in range(n):
            m.readout[w] = (0.05, 0.06)
    elif scenario == "combined":     # everything moderate
        for w in range(n):
            m.t1_us[w] = 30.0
            m.t2_us[w] = 18.0
            m.readout[w] = (0.02, 0.025)
        m.gate_infidelity = {"1q": 0.001, "cx": 0.006, "cz": 0.006,
                             "swap": 0.018, "cp": 0.006}
        drift = 0.0002
    else:
        raise ValueError(scenario)
    for w in range(n):
        m.drift_rate[w] = drift
    return m


SCENARIO_EPS = {"coherent": 0.15, "decoherence": 0.01,
                "readout": 0.01, "combined": 0.05}


def main() -> int:
    K_twirl = 8
    all_rows = []
    failures = []
    t0 = time.perf_counter()
    for scenario in ("coherent", "decoherence", "readout", "combined"):
        eps = SCENARIO_EPS[scenario]
        print(f"=== scenario: {scenario} (eps_coh={eps}) ===")
        # families with meaningful deterministic-outcome metrics; the QFT
        # family has a (near-)uniform ideal distribution under this harness
        # and is excluded rather than reported with a vacuous metric
        for family in ("ghz", "bv"):
            for n in (3, 4, 5):
                circ = build(family, n)
                nz = scenario_model(n, scenario)
                mask = good_mask_for(family, n, circ=circ)

                def run_pipeline(c, mem=False, twirl=False):
                    # every pipeline is averaged over the SAME K noise
                    # instances; twirled pipelines pair instance k with
                    # twirl variant k and place DD on THAT variant's own
                    # schedule (DD last - it must see the final idle
                    # structure, and Pauli-merging must never eat it)
                    acc = None
                    for k in range(K_twirl):
                        if twirl:
                            cc = insert_dd(pauli_twirl(c, seed=k), nz,
                                           min_idle_ns=100)
                        else:
                            cc = c
                        rho = run_noisy(cc, nz, eps, seed=k)
                        d = np.real(np.diag(rho))
                        acc = d if acc is None else acc + d
                    acc = np.clip(acc / K_twirl, 0, None)
                    acc /= acc.sum()
                    return success(np.diag(acc).astype(complex), mask,
                                   nz, mem=mem)

                rows = {}
                rows["raw"] = run_pipeline(circ)
                opt = optimize_search(circ, verify=False)
                opt = expand_for_suppression(opt)
                rows["optimized"] = run_pipeline(opt)
                rows["+dd"] = run_pipeline(insert_dd(opt, nz, min_idle_ns=100))
                rows["+dd+twirl"] = run_pipeline(opt, twirl=True)
                rows["+dd+twirl+mem"] = run_pipeline(opt, twirl=True, mem=True)
                factor = rows["+dd+twirl+mem"] / max(rows["raw"], 1e-12)
                all_rows.append({"scenario": scenario, "family": family,
                                 "n": n, **rows, "factor": round(factor, 2)})
                cells = " ".join(f"{k}={v:.4f}" for k, v in rows.items())
                verdict = "PASS" if factor > 1.0 else "FAIL"
                if factor <= 1.0:
                    failures.append((scenario, family, n, factor))
                print(f"  {family:4s} n={n}: {cells} | x{factor:6.2f} {verdict}")
        print()
    dt = time.perf_counter() - t0
    doc = {"schema": "qproof-suppress/3", "twirl_variants": K_twirl,
           "gate": "full pipeline factor > 1.0 in every scenario",
           "wall_time_s": round(dt, 1), "results": all_rows}
    (REPO / "suppress_results.json").write_text(json.dumps(doc, indent=2) + chr(10),
                                                encoding="utf-8")
    print(f"wrote suppress_results.json ({dt:.0f}s)")
    if failures:
        print(f"GATE FAILED: {len(failures)} scenario(s) at or below raw:")
        for s, f, n, x in failures:
            print(f"  {s}/{f}/n={n}: x{x}")
        return 1
    print("GATE PASSED: full pipeline beats raw in every scenario.")
    return 0


def dd_sequence_comparison() -> int:
    """DD sequence-family comparison (suppression depth phase).

    For the two DD-relevant scenarios, the full suppression pipeline is
    run once per decoupling sequence (xy4 / xy8 / xzx / pdd4 / auto).
    Gate: every sequence's pipeline beats raw in every cell of the
    DD-relevant scenarios - sequence choice may differ in degree but
    none may lose."""
    K = 8
    all_rows = []
    failures = []
    for scenario in ("decoherence", "combined"):
        eps = SCENARIO_EPS[scenario]
        print(f"=== dd-sequences under {scenario} ===")
        for family in ("ghz", "bv"):
            for n in (4, 5):
                circ = build(family, n)
                nz = scenario_model(n, scenario)
                mask = good_mask_for(family, n, circ=circ)
                opt = optimize_search(circ)
                opt = expand_for_suppression(opt)

                def run_seq(seq):
                    acc = None
                    for k in range(K):
                        cc = insert_dd(pauli_twirl(opt, seed=k), nz,
                                       min_idle_ns=100, sequence=seq)
                        rho = run_noisy(cc, nz, eps, seed=k)
                        d = np.real(np.diag(rho))
                        acc = d if acc is None else acc + d
                    acc = np.clip(acc / K, 0, None)
                    acc /= acc.sum()
                    return success(np.diag(acc).astype(complex), mask,
                                   nz, mem=True)

                raw_acc = None
                for k in range(K):
                    rho = run_noisy(circ, nz, eps, seed=k)
                    d = np.real(np.diag(rho))
                    raw_acc = d if raw_acc is None else raw_acc + d
                raw_acc = np.clip(raw_acc / K, 0, None)
                raw_acc /= raw_acc.sum()
                raw_p = success(np.diag(raw_acc).astype(complex), mask,
                                nz, mem=False)

                seq_row = {}
                for seq in ("xy4", "xy8", "xzx", "pdd4", "auto"):
                    seq_row[seq] = run_seq(seq)
                worst = min(seq_row.values())
                row = {"scenario": scenario, "family": family, "n": n,
                       "raw": round(raw_p, 4),
                       **{s: round(v, 4) for s, v in seq_row.items()}}
                all_rows.append(row)
                verdict = "PASS" if worst > raw_p else "FAIL"
                if verdict == "FAIL":
                    failures.append(row)
                cells = " ".join(f"{s}={v:.4f}" for s, v in seq_row.items())
                print(f"  {family:4s} n={n}: raw={raw_p:.4f} {cells} {verdict}",
                      flush=True)
    (Path(__file__).resolve().parents[1] / "dd_sequence_results.json").write_text(
        json.dumps({"schema": "qproof-ddseq/1", "results": all_rows},
                   indent=2) + "\n", encoding="utf-8")
    if failures:
        print(f"SEQUENCE GATE FAILED: {failures}")
        return 1
    print("SEQUENCE GATE PASSED: every DD sequence beats raw everywhere.")
    return 0


if __name__ == "__main__":
    if "--dd-sequences" in sys.argv:
        sys.exit(dd_sequence_comparison())
    sys.exit(main())
