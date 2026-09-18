"""Head-to-head: Compact Suppression vs Qiskit-O3 vs raw execution.

The protocol mirrors how closed suppression services market themselves
("our pipeline vs default Qiskit execution") but open and reproducible:

  raw            : the input circuit, executed as-is
  qiskit-O3      : qiskit transpile(optimization_level=3), same basis
  qiskit-O3+supp : the qiskit-O3 circuit + the FULL Compact suppression
                   stack (twirl x K, benefit-gated DD, readout mitigation)
  compact-full   : compactq optimize + the same full suppression stack

Every row runs through the SAME density-matrix noise simulator with the
same seeds, coherent axes and readout confusion; the metric is the
probability of the ideal outcome.  Suppression gates (benefit-gated DD,
exactness-proven twirl) compose with EITHER compiler - the point of the
qiskit-O3+supp row.

Usage: python scripts/head_to_head.py   (writes head_to_head_results.json)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


def git_sha():
    import subprocess
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, cwd=str(REPO)).stdout.strip()[:12]
    except Exception:
        return "unknown"


from datetime import datetime, timezone  # noqa: E402
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from compactq import Circuit, Gate, optimize_search  # noqa: E402
from compactq.noise import NoiseModel  # noqa: E402
from compactq.suppress import (pauli_twirl, insert_dd,  # noqa: E402
                               expand_for_suppression)
from compactq.mitigate import mitigate_counts  # noqa: E402

from suppress_bench import (run_noisy, success, good_mask_for,  # noqa: E402
                            scenario_model, SCENARIO_EPS)

K = 8


def qiskit_o3(circ: Circuit) -> Circuit:
    """qiskit optimization_level=3 on the same qubits, translated back."""
    from qiskit import transpile
    from compactq.qiskit_bridge import to_qiskit, from_qiskit
    qc = to_qiskit(circ)
    out = transpile(qc, optimization_level=3,
                    basis_gates=["rz", "sx", "x", "h", "cx", "cz"],
                    coupling_map=None, seed_transpiler=7)
    back = from_qiskit(out)
    assert back.num_qubits == circ.num_qubits
    return back


def full_stack(circ: Circuit, nz: NoiseModel, eps: float, mask):
    """optimize -> expand -> [twirl variant k -> DD] per variant -> MEM,
    returning a callable-free pipeline evaluated by the caller's runner."""
    opt = optimize_search(circ, verify=False)
    opt = expand_for_suppression(opt)

    def run(mem: bool = False) -> float:
        acc = None
        for k in range(K):
            cc = insert_dd(pauli_twirl(opt, seed=k), nz, min_idle_ns=100)
            rho = run_noisy(cc, nz, eps, seed=k)
            d = np.real(np.diag(rho))
            acc = d if acc is None else acc + d
        acc = np.clip(acc / K, 0, None)
        acc /= acc.sum()
        return success(np.diag(acc).astype(complex), mask, nz, mem=mem)
    return run


def plain(circ: Circuit, nz: NoiseModel, eps: float, mask) -> float:
    acc = None
    for k in range(K):
        rho = run_noisy(circ, nz, eps, seed=k)
        d = np.real(np.diag(rho))
        acc = d if acc is None else acc + d
    acc = np.clip(acc / K, 0, None)
    acc /= acc.sum()
    return success(np.diag(acc).astype(complex), mask, nz, mem=False)


def qaoa_ring(n: int) -> Circuit:
    """One QAOA layer on a ring: cost rzz(pi/4) expanded CX-style, full
    rx mixer - a parameterized hybrid subroutine like Fire Opal's solvers."""
    ops = [Gate("h", (), (w,)) for w in range(n)]
    for w in range(n):
        a, b = w, (w + 1) % n
        ops += [Gate("cx", (), (a, b)), Gate("rz", (np.pi / 2,), (b,)),
                Gate("cx", (), (a, b))]
    ops += [Gate("rx", (np.pi / 3,), (w,)) for w in range(n)]
    return Circuit(n, ops)


def build(family: str, n: int) -> Circuit:
    if family == "ghz":
        ops = [Gate("h", (), (0,))]
        ops += [Gate("cx", (), (j, j + 1)) for j in range(n - 1)]
        return Circuit(n, ops)
    if family == "bv":
        ops = [Gate("x", (), (n - 1,))]
        ops += [Gate("h", (), (j,)) for j in range(n - 1)]
        ops += [Gate("cz", (), (j, n - 1)) for j in range(n - 1)]
        ops += [Gate("h", (), (j,)) for j in range(n - 1)]
        return Circuit(n, ops)
    if family == "qaoa":
        return qaoa_ring(n)
    raise ValueError(family)


def main() -> int:
    t0 = time.perf_counter()
    rows = []
    losses = []
    sums = {}   # scenario -> [sum raw, sum qiskit-o3, sum compact]
    for scenario in ("coherent", "decoherence", "combined"):
        nz_factory = scenario_model
        eps = SCENARIO_EPS[scenario]
        print(f"=== {scenario} (eps_coh={eps}) ===")
        for family in ("ghz", "bv", "qaoa"):
            for n in (4, 6):
                circ = build(family, n)
                nz = nz_factory(n, scenario)
                mask = good_mask_for(family, n, circ=circ)

                row = {"scenario": scenario, "family": family, "n": n}
                row["raw"] = plain(circ, nz, eps, mask)

                q3 = qiskit_o3(circ)
                row["qiskit-o3"] = plain(q3, nz, eps, mask)
                # qiskit-O3 + the suppression stack (portability check)
                acc = None
                for k in range(K):
                    cc = insert_dd(pauli_twirl(q3, seed=k), nz, min_idle_ns=100)
                    rho = run_noisy(cc, nz, eps, seed=k)
                    d = np.real(np.diag(rho))
                    acc = d if acc is None else acc + d
                acc = np.clip(acc / K, 0, None)
                acc /= acc.sum()
                row["qiskit-o3+supp"] = success(np.diag(acc).astype(complex),
                                                mask, nz, mem=True)

                run_full = full_stack(circ, nz, eps, mask)
                row["compact-full"] = run_full(mem=True)

                # gate: suppression must pay for itself vs RAW in every
                # cell, and the pipeline must beat the qiskit-O3 compiler
                # baseline in aggregate per scenario.  (Single cells can
                # legitimately favor bare compilation: RC trades coherent
                # error cancellation for stochastic robustness - visible
                # for BOTH stacks as qiskit+supp < qiskit-O3.)
                if row["compact-full"] <= row["raw"]:
                    losses.append((scenario, family, n,
                                   row["compact-full"], row["raw"]))
                sums.setdefault(scenario, [0.0, 0.0, 0.0])
                sums[scenario][0] += row["raw"]
                sums[scenario][1] += row["qiskit-o3"]
                sums[scenario][2] += row["compact-full"]
                rows.append(row)
                print(f"  {family:5s} n={n}: raw={row['raw']:.4f} "
                      f"qiskit-O3={row['qiskit-o3']:.4f} "
                      f"qiskit+supp={row['qiskit-o3+supp']:.4f} "
                      f"compact={row['compact-full']:.4f}")
    dt = time.perf_counter() - t0
    print("\nscenario aggregates (mean success probability):")
    for s, (r, q, c) in sums.items():
        cnt = sum(1 for row in rows if row["scenario"] == s)
        print(f"  {s:12s} raw={r/cnt:.4f} qiskit-O3={q/cnt:.4f} "
              f"compact-full={c/cnt:.4f}")
    doc = {"schema": "qproof-headtohead/1", "twirl_variants": K,
           "generated": datetime.now(timezone.utc).isoformat(),
           "commit": git_sha(),
           "wall_time_s": round(dt, 1), "results": rows,
           "aggregates": {s: {"raw": r / sum(1 for row in rows if row["scenario"] == s),
                              "qiskit-o3": q / sum(1 for row in rows if row["scenario"] == s),
                              "compact-full": c / sum(1 for row in rows if row["scenario"] == s)}
                          for s, (r, q, c) in sums.items()}}
    (REPO / "head_to_head_results.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    if losses:
        print(f"GATE FAILED: compact-full at or below raw in {len(losses)} cell(s):")
        for s, f, n, c, b in losses:
            print(f"  {s}/{f}/n={n}: compact={c:.4f} raw={b:.4f}")
        return 1
    for s, (r, q, c) in sums.items():
        if c <= q:
            print(f"GATE FAILED: compact-full does not beat qiskit-O3 "
                  f"in aggregate for {s}")
            return 1
    print("GATE PASSED: compact-full > raw in every cell, "
          "> qiskit-O3 in aggregate in every scenario.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
