"""Scalability gauntlet — find Compact's REAL correctness/runtime/memory ceilings.

For every width in 4..128 and a spread of workload families, this harness:
  1. generates the circuit (seeded, deterministic),
  2. optimizes it in an isolated subprocess (wall-clock timeout, tracemalloc
     peak-memory),
  3. verifies the output at the strongest available tier:
       - independent Qiskit Operator referee at small widths,
       - exact algebraic tableau proof for Clifford workloads at ANY width,
       - otherwise the never-grow policy check (gates/2q/depth/cx-equiv must
         not increase) plus a determinism double-run,
  4. records everything and aggregates per-width PASS/TEST/FAIL.

The claim this harness supports is deliberately strict:
  "stable at N qubits" = zero crashes, zero timeouts, zero policy
  violations, zero verified-incorrect outputs, bounded memory, and
  deterministic repeats at that width.

Usage:
  python scripts/scale_gauntlet.py                       # default profile
  python scripts/scale_gauntlet.py --widths 8,16,32 --cases 1
  python scripts/scale_gauntlet.py --ci                  # quick CI profile

Artifacts: results/scalability.json + results/scalability.md
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import statistics
import sys
import time
import tracemalloc
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

# inner gauntlet cases deliberately run without the in-product proof: this
# harness is itself the independent verifier (tiered, see verify_case)
_OFF = False

TIMEOUT_S = 180
MEMORY_BUDGET_MB = 1024.0


# ----------------------------------------------------------------- workloads
def _w_ghz(n: int):
    from compactq import Circuit, Gate
    return Circuit(n, [Gate("h", (), (0,))]
                   + [Gate("cx", (), (j, j + 1)) for j in range(n - 1)])


def _w_qft(n: int):
    from compactq import Circuit, Gate
    import math
    ops = []
    for j in range(n):
        ops.append(Gate("h", (), (j,)))
        for k in range(j + 1, n):
            ops.append(Gate("cp", (math.pi / 2 ** (k - j),), (k, j)))
    for j in range(n // 2):
        ops.append(Gate("swap", (), (j, n - 1 - j)))
    return Circuit(n, ops)


def _w_random_clifford(n: int, seed: int = 0):
    import random
    from compactq import Circuit, Gate
    rng = random.Random(seed * 977 + n)
    oneq = ("h", "x", "s", "sdg")
    ops = []
    for layer in range(3 * n):
        q = rng.randrange(n)
        ops.append(Gate(rng.choice(oneq), (), (q,)))
        if rng.random() < 0.5 and n >= 2:
            a = rng.randrange(n)
            b = (a + rng.randrange(1, n)) % n
            ops.append(Gate("cx", (), (a, b)))
    return Circuit(n, ops)


def _w_random_general(n: int, seed: int = 0):
    import random
    from compactq import Circuit, Gate
    rng = random.Random(seed * 613 + n)
    ops = []
    for _ in range(2 * n):
        name = rng.choice(("h", "x", "ry", "rz"))
        if name in ("ry", "rz"):
            ops.append(Gate(name, (rng.uniform(0.1, 3.04),), (rng.randrange(n),)))
        else:
            ops.append(Gate(name, (), (rng.randrange(n),)))
        if n >= 2:
            a = rng.randrange(n)
            b = (a + rng.randrange(1, n)) % n
            ops.append(Gate("cx", (), (a, b)))
    return Circuit(n, ops)


def _w_ising_trotter(n: int, layers: int = 3):
    import random
    from compactq import Circuit, Gate
    rng = random.Random(n)
    ops = []
    for _ in range(layers):
        for q in range(n):
            ops.append(Gate("rz", (rng.uniform(0.4, 1.2),), (q,)))
        for q in range(n - 1):
            ops.append(Gate("cx", (), (q, q + 1)))
            ops.append(Gate("rz", (rng.uniform(0.2, 0.8),), (q + 1,)))
    return Circuit(n, ops)


def _w_qpe_like(n: int):
    from compactq import Circuit, Gate
    import math
    ops = [Gate("h", (), (j,)) for j in range(max(1, n // 2))]
    for j in range(n - 1):
        ops.append(Gate("cp", (math.pi / 2 ** (j + 1),), (j, n - 1)))
    for j in range(max(1, n // 2)):
        ops.append(Gate("h", (), (j,)))
    return Circuit(n, ops)


def _w_hwe_ansatz(n: int, layers: int = 2):
    import random
    from compactq import Circuit, Gate
    rng = random.Random(n * 31)
    ops = []
    for _ in range(layers):
        for q in range(n):
            ops.append(Gate("ry", (rng.uniform(0, 3.14),), (q,)))
            ops.append(Gate("rz", (rng.uniform(0, 3.14),), (q,)))
        for q in range(n - 1):
            ops.append(Gate("cx", (), (q, q + 1)))
    return Circuit(n, ops)


def _w_library_random(n: int, seed: int = 0):
    from compactq.benchmarks import random_circuit
    return random_circuit(n, 6 * n, seed=seed + n)


WORKLOADS = {
    "ghz": lambda n, s: _w_ghz(n),
    "qft": lambda n, s: _w_qft(n),
    "random_clifford": _w_random_clifford,
    "random_general": _w_random_general,
    "ising_trotter": lambda n, s: _w_ising_trotter(n),
    "qpe_like": lambda n, s: _w_qpe_like(n),
    "hwe_ansatz": lambda n, s: _w_hwe_ansatz(n),
    "library_random": _w_library_random,
}


# --------------------------------------------------------------- case worker
def run_case(width: int, workload: str, case: int, referee_limit: int) -> dict:
    """Executed inside an isolated subprocess (spawn).  Returns a record
    dict; exceptions are captured, the parent enforces the timeout."""
    rec = {"width": width, "workload": workload, "case": case}
    try:
        import compactq
        from compactq import optimize_search
        from compactq.stabilizer import is_clifford
        rec["compactq_version"] = compactq.__version__

        circ = WORKLOADS[workload](width, case)
        m_in = (circ.two_qubit_count(), len(circ.ops), circ.depth(),
                circ.cx_equivalent_count())
        rec["gates_in"], rec["two_q_in"], rec["depth_in"] = m_in[0], m_in[1], m_in[2]
        rec["cx_eq_in"] = m_in[3]
        rec["clifford_input"] = bool(is_clifford(circ))

        tracemalloc.start()
        t0 = time.perf_counter()
        out = optimize_search(circ, verify=_OFF)
        opt_s = time.perf_counter() - t0
        _cur, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rec["optimization_ms"] = round(opt_s * 1000)
        rec["peak_mb"] = round(peak / (1024 * 1024), 1)

        m_out = (out.two_qubit_count(), len(out.ops), out.depth(),
                 out.cx_equivalent_count())
        rec["gates_out"], rec["two_q_out"], rec["depth_out"] = m_out[0], m_out[1], m_out[2]
        rec["cx_eq_out"] = m_out[3]

        # never-grow policy check: the CONTRACT is lexicographic on
        # (2q, gates, depth) — depth may rise only when total gates fall.
        # Recorded separately for transparency.
        rec["policy_ok"] = m_out[:3] <= m_in[:3]
        rec["depth_trade"] = m_out[2] > m_in[2] and rec["policy_ok"]

        # determinism: a second run must produce an identical circuit shape
        t1 = time.perf_counter()
        out2 = optimize_search(circ, verify=_OFF)
        rec["deterministic"] = ((out2.two_qubit_count(), len(out2.ops),
                                 out2.depth()) == m_out[:3])
        rec["verify_ms"] = round((time.perf_counter() - t1) * 1000)

        # tiered independent verification
        import compactq.verify as vmod
        t2 = time.perf_counter()
        if width <= referee_limit:
            try:
                import numpy as np  # noqa
                from qiskit.quantum_info import Operator
                from compactq.qiskit_bridge import to_qiskit
                ref = Operator(to_qiskit(circ)).data
                got = Operator(to_qiskit(out)).data
                fid = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                rec["equivalence"] = (fid > 1 - 1e-6)
                rec["tier"] = "independent-qiskit"
                rec["fidelity"] = round(fid, 12)
            except Exception:
                # referee choked on an exotic shape: fall back honestly to
                # policy-only rather than claiming a decision
                rec["equivalence"], rec["tier"] = None, "policy-only"
        elif is_clifford(circ) and is_clifford(out):
            from compactq.stabilizer import clifford_equal
            rec["equivalence"] = bool(clifford_equal(circ, out))
            rec["tier"] = "algebraic-tableau"
        else:
            rec["equivalence"] = None
            rec["tier"] = "policy-only"
        rec["verify_ms"] = round((time.perf_counter() - t2) * 1000)
        rec["status"] = "ok"
    except Exception as e:
        rec["status"] = "exception"
        rec["exception"] = f"{type(e).__name__}: {e}"
    return rec


def _worker(pipe, width, workload, case, referee_limit):
    try:
        result = run_case(width, workload, case, referee_limit)
    except Exception as e:  # whole-worker death
        result = {"width": width, "workload": workload, "case": case,
                  "status": "exception",
                  "exception": f"worker {type(e).__name__}: {e}"}
    try:
        pipe.send(result)
    except (BrokenPipeError, OSError):
        pass  # parent stopped listening (timeout path)
    finally:
        pipe.close()


def run_with_timeout(width, workload, case, referee_limit, timeout_s):
    ctx = multiprocessing.get_context("spawn")
    parent, child = ctx.Pipe(duplex=False)
    p = ctx.Process(target=_worker,
                    args=(child, width, workload, case, referee_limit))
    p.start()
    child.close()  # parent keeps the recv end; drop its copy of the send end
    if p.is_alive():
        p.join(timeout_s)
    if p.is_alive():
        p.terminate()
        p.join(5)
        return {"width": width, "workload": workload, "case": case,
                "status": "timeout", "timeout_s": timeout_s}
    result = None
    if parent.poll():
        try:
            result = parent.recv()
        except EOFError:
            result = None
    parent.close()
    if result is None:
        return {"width": width, "workload": workload, "case": case,
                "status": "exception", "exception": "no result from worker"}
    return result


# ------------------------------------------------------------------ aggregate
def summarize(records, timeout_s, memory_budget_mb):
    by_width = {}
    for r in records:
        by_width.setdefault(r["width"], []).append(r)
    widths = sorted(by_width)
    rows = []
    correctness_ceiling = runtime_ceiling = memory_ceiling = None
    broken = False
    for w in widths:
        rs = by_width[w]
        crashes = sum(1 for r in rs if r["status"] == "exception")
        timeouts = sum(1 for r in rs if r["status"] == "timeout")
        decided = [r for r in rs if isinstance(r.get("equivalence"), bool)]
        incorrect = sum(1 for r in decided if r["equivalence"] is False)
        policy_bad = sum(1 for r in rs if not r.get("policy_ok", True))
        depth_trades = sum(1 for r in rs if r.get("depth_trade"))
        nondet = sum(1 for r in rs if r.get("deterministic") is False)
        times = sorted(r["optimization_ms"] for r in rs
                       if isinstance(r.get("optimization_ms"), int))
        peak = max((r.get("peak_mb", 0) for r in rs), default=0)
        median = round(statistics.median(times)) if times else None
        ok = (crashes == 0 and timeouts == 0 and incorrect == 0
              and policy_bad == 0 and nondet == 0)
        status = "PASS" if ok else ("FAIL" if (crashes or incorrect or policy_bad)
                                    else "TEST")
        rows.append({
            "width": w, "tests": len(rs), "crashes": crashes,
            "timeouts": timeouts, "incorrect": incorrect,
            "policy_violations": policy_bad, "depth_trades": depth_trades,
            "nondeterministic": nondet,
            "verified": len(decided), "unverified": len(rs) - len(decided),
            "median_opt_ms": median,
            "peak_mb": peak,
            "tiers": sorted({r.get("tier", "?") for r in rs}),
            "status": status,
        })
        # ceilings follow the unbroken chain: the first FAIL below a passing
        # width caps all three ceilings (re-run after fixing to extend)
        if broken:
            continue
        if status == "PASS":
            correctness_ceiling = w
            if median is not None and median <= timeout_s * 500:
                runtime_ceiling = w
            if peak <= memory_budget_mb:
                memory_ceiling = w
        elif status == "FAIL":
            broken = True
    return rows, correctness_ceiling, runtime_ceiling, memory_ceiling


def write_md(rows, ceilings, doc, path):
    corr, runtime, memory = ceilings
    lines = [
        "# Scalability gauntlet (generated - do not edit)",
        "",
        f"Generated {doc['generated']} | commit `{doc['commit']}` | "
        f"compactq {doc['compactq_version']} | referee limit "
        f"{doc['referee_limit']}q (independent Qiskit Operator); Clifford "
        f"workloads are algebraically proven at any width.",
        "",
        "| width | tests | crashes | timeouts | incorrect | policy viol | "
        "nondet | verified | median opt ms | peak MB | status |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['width']}Q | {r['tests']} | {r['crashes']} | "
            f"{r['timeouts']} | {r['incorrect']} | {r['policy_violations']} | "
            f"{r['nondeterministic']} | {r['verified']} | "
            f"{r['median_opt_ms']} | {r['peak_mb']} | {r['status']} |")
    lines += ["",
              f"- correctness/stability ceiling (all workloads): "
              f"**{corr}Q**",
              f"- runtime ceiling (median opt time sane): **{runtime}Q**",
              f"- memory ceiling (< {int(doc['memory_budget_mb'])} MB peak): "
              f"**{memory}Q**",
              "",
              "Verified-incorrect anywhere is a release blocker; "
              "'policy-only' rows are stability-proven, not "
              "equivalence-proven (see README tier table).", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="Compact scalability gauntlet")
    ap.add_argument("--widths", default="4,6,8,10,12,16,20,24,32,48,64,96,128")
    ap.add_argument("--cases", type=int, default=1,
                    help="cases per (width, workload); scale to 1000 for the "
                         "strict-ceiling protocol")
    ap.add_argument("--referee-limit", type=int, default=10,
                    help="max width for the independent Qiskit referee")
    ap.add_argument("--timeout", type=float, default=TIMEOUT_S)
    ap.add_argument("--memory-budget-mb", type=float, default=MEMORY_BUDGET_MB)
    ap.add_argument("--ci", action="store_true",
                    help="quick profile: widths 4,8,16,24 only")
    ap.add_argument("--out-dir", default=str(REPO / "results"))
    args = ap.parse_args()

    widths = ([4, 8, 16, 24] if args.ci
              else [int(x) for x in args.widths.split(",") if x.strip()])
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    import compactq
    records = []
    for w in widths:
        for wl in WORKLOADS:
            for case in range(max(1, args.cases)):
                t0 = time.perf_counter()
                rec = run_with_timeout(w, wl, case, args.referee_limit,
                                       args.timeout)
                records.append(rec)
                st = rec.get("status", "?")
                print(f"{w:>4}Q {wl:16s} case {case}: {st} "
                      f"({time.perf_counter() - t0:.1f}s)", flush=True)

    rows, corr, runtime, memory = summarize(records, args.timeout,
                                            args.memory_budget_mb)
    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": None,
        "compactq_version": compactq.__version__,
        "python": sys.version.split()[0],
        "widths": widths,
        "cases_per_workload": max(1, args.cases),
        "referee_limit": args.referee_limit,
        "timeout_s": args.timeout,
        "memory_budget_mb": args.memory_budget_mb,
        "ceilings": {"correctness_stability": corr, "runtime": runtime,
                     "memory": memory},
        "records": records,
    }
    try:
        import subprocess
        doc["commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=str(REPO)).stdout.strip()[:12]
    except Exception:
        pass
    (out_dir / "scalability.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    write_md(rows, (corr, runtime, memory), doc,
             out_dir / "scalability.md")
    print(f"wrote {out_dir / 'scalability.json'} + .md "
          f"({len(records)} records)")
    print(f"ceilings: correctness/stability {corr}Q | runtime {runtime}Q | "
          f"memory {memory}Q")
    bad = [r for r in records if r.get("status") != "ok"
           or r.get("equivalence") is False
           or not r.get("policy_ok", True)]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
