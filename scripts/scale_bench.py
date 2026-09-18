"""Scale benchmark: suppression layers at n = 16..24 (stabilizer regime).

The density-matrix study (scripts/suppress_bench.py) tops out near 6
qubits; this harness pushes the same pipeline to the scale regime where
closed suppression services show their headline demos - using the CHP
stabilizer simulator (compactq.stabsim).  Honest scope: stochastic gate
noise + readout confusion (coherent rotations and amplitude damping are
not Pauli-representable; the small-n studies cover those).

Device: a line with ONE POISONED edge near its end (20% cx error vs 1%
elsewhere).  Compared mappings:
  naive-map : identity placement, blind routing (no calibration sight)
  aware-map : fidelity-aware placement + poison-sighted routing
  aware+mem : aware-map + tensored readout inversion

Gate: the full stack (aware+mem) beats the naive mapping in every cell.
Reproducible, seed-fixed.
"""
from __future__ import annotations

import json
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from compactq import Circuit, Gate, mitigate_counts, stab_sample  # noqa: E402
from compactq.hardware import layout_aware, route_aware  # noqa: E402
from compactq.noise import NoiseModel  # noqa: E402
from compactq.stabsim import stab_run_shot  # noqa: E402
from compactq.target import Target  # noqa: E402

SHOTS = 2000
P_CX = 0.01          # uniform two-qubit depolarizing rate
P_1Q = 0.002         # one-qubit depolarizing rate
RO = 0.03            # per-wire readout confusion (symmetric)
POISON_FID = 0.80    # the poisoned edge's cx fidelity (vs 0.99)


def sample_with_readout(circ, shots, rng, ro_map, noise=None):
    out = {}
    for _ in range(shots):
        s = list(stab_run_shot(circ, rng, noise=noise))
        for w in range(circ.num_qubits):
            if rng.random() < ro_map[w]:
                s[w] = '1' if s[w] == '0' else '0'
        key = ''.join(s)
        out[key] = out.get(key, 0) + 1
    return out


def ghz(n):
    return Circuit(n, [Gate('h', (), (0,))]
                   + [Gate('cx', (), (j, j + 1)) for j in range(n - 1)])


def bv(n):
    return Circuit(n, [Gate('x', (), (n - 1,))]
                   + [Gate('h', (), (j,)) for j in range(n - 1)]
                   + [Gate('cz', (), (j, n - 1)) for j in range(n - 1)]
                   + [Gate('h', (), (j,)) for j in range(n - 1)])


def ideal_mass(counts, support):
    tot = sum(counts.values()) or 1
    return sum(c for k, c in counts.items() if k in support) / tot


def main() -> int:
    t0 = time.perf_counter()
    rows = []
    failures = []
    for family, build in (("ghz", ghz), ("bv", bv)):
        for n in (16, 24):
            device = n + 4
            coupling = [(j, j + 1) for j in range(device - 1)]
            poisoned = (1, 2)   # near the end: avoidable by placement

            circ = build(n)
            support = {'0' * n, '1' * n} if family == 'ghz' else {'1' * n}
            rng = random.Random(100 + n)

            # calibration sight for mapping/routing
            noise = NoiseModel(device)
            noise.gate_infidelity = {'cx': P_CX, '1q': P_1Q,
                                     ('cx', poisoned): 1.0 - POISON_FID}
            target = Target(cx_fidelity={frozenset(poisoned): POISON_FID},
                            default_cx_fidelity=1.0 - P_CX)

            def map_and_sample(placement, sighted):
                dev = Circuit(device, [Gate(g.name, g.params,
                                            tuple(placement[w]
                                                  for w in g.qubits))
                                      for g in circ.ops])
                routed, final_pos = route_aware(
                    dev, coupling, target=(target if sighted else None),
                    restore=False)
                ro_map = {w: RO for w in range(device)}
                counts = sample_with_readout(routed, SHOTS, rng, ro_map,
                                             noise=noise)
                logical = {}
                for k, c in counts.items():
                    lk = ''.join(k[final_pos[placement[w]]]
                                 for w in range(n))
                    logical[lk] = logical.get(lk, 0) + c
                return ideal_mass(logical, support), logical

            identity = {w: w for w in range(n)}
            naive_p, _ = map_and_sample(identity, sighted=False)
            _, placement = layout_aware(circ, coupling, noise,
                                        max_placements=200000)
            aware_p, aware_counts = map_and_sample(placement, sighted=True)
            mit = mitigate_counts(aware_counts,
                                  {w: (RO, RO) for w in range(n)}, n)
            mem_p = sum(v for k, v in mit.items() if k in support)

            row = {"family": family, "n": n,
                   "naive-map": round(naive_p, 4),
                   "aware-map": round(aware_p, 4),
                   "aware+mem": round(mem_p, 4)}
            rows.append(row)
            verdict = "PASS" if mem_p > naive_p and aware_p > naive_p - 0.05 \
                else "FAIL"
            if verdict == "FAIL":
                failures.append(row)
            print(f"{family:4s} n={n}: naive={naive_p:.4f} "
                  f"aware={aware_p:.4f} aware+mem={mem_p:.4f} {verdict}",
                  flush=True)

    from provenance import environment
    doc = {"schema": "qproof-scale/1", "shots": SHOTS,
           "p_cx": P_CX, "readout": RO,
           "generated": datetime.now(timezone.utc).isoformat(),
           "commit": __import__("provenance").git_sha(),
           "environment": environment(),
           "wall_time_s": round(time.perf_counter() - t0, 1),
           "results": rows}
    (REPO / "scale_results.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    if failures:
        print(f"GATE FAILED: {failures}")
        return 1
    print("GATE PASSED: aware+mem beats the naive mapping at every scale.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
