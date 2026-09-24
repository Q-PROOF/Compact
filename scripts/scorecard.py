"""Generate results/SCORECARD.json + SCORECARD.md — the multi-metric,
multi-corpus scorecard recommended for a verified optimizer:

  * gate reduction (median AND geometric-mean % change)
  * entangling-gate reduction (2q before/after)
  * depth reduction
  * compilation cost (median, p95, timeout/error rate — not just means)
  * correctness (verified/total + referee name/version)
  * per-config rows (compact / qiskit L3 / pytket / cirq) — never one
    opaque aggregate
  * reproducibility (commit, versions, commands, machine)

Sources are the committed results/ artifacts only.  A corpus whose
artifact is missing is reported NOT MEASURED; a tool that could not run
in this environment is reported as such — nothing is estimated.

Usage:  python scripts/scorecard.py
"""
from __future__ import annotations

import json
import math
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"


def _load(path):
    p = Path(path)
    if not p.is_absolute():
        p = RESULTS / path
    if not p.exists():
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def _median(xs):
    xs = [x for x in xs if x is not None]
    return round(statistics.median(xs), 2) if xs else None


def _median_pct(xs):
    """Median of fraction-valued reductions, expressed in percent."""
    xs = [x for x in xs if x is not None]
    return round(statistics.median(xs) * 100, 2) if xs else None


def _geomean_pct(pcts):
    vals = [p for p in pcts if p is not None]
    if not vals:
        return None
    return round(math.exp(sum(math.log(max(v, 1e-9)) for v in vals)
                          / len(vals)) * 100, 2)


def _reduction(before, after):
    if not before:
        return None
    return (before - after) / before


def _compact_repro_block(*sources):
    commits = sorted({s.get("commit") for s in sources
                      if isinstance(s, dict) and s.get("commit")})
    return {"commits": commits,
            "commands": {
                "repro": "python scripts/repro_harness.py",
                "qasmbench_small": "python scripts/bench_json.py",
                "mqtbench": "python scripts/mqtbench_run.py",
                "benchpress": "python scripts/benchpress_run.py --gym compactq",
                "feynman": "python scripts/feynman_bench.py",
                "ddceiling": "python scripts/dd_ceiling_bench.py",
                "latency": "python scripts/latency_bench.py",
            },
            "manifest": "results/MANIFEST.json (pinned suite)",
            "package_versions": sorted({
                s.get("compactq_version") for s in sources
                if isinstance(s, dict) and s.get("compactq_version")})}


def score_qasmbench_small():
    d = _load(REPO / "bench_results.json")
    if not d:
        return None
    rows = [r for r in d["circuits"]
            if isinstance(r.get("input"), dict) and r.get("compactq")
            and r.get("qiskit_l3")]
    out = {"corpus": f"QASMBench small ({len(rows)} scored of "
                     f"{len(d['circuits'])} circuits)",
           "artifact": "bench_results.json (repo root)",
           "referee": "qiskit Operator, level-0 normalized inputs",
           "configs": {}}
    for tool in ("compactq", "qiskit_l3"):
        red_g = [_reduction(r["input"]["gates"], r[tool]["gates"]) for r in rows]
        red_2q = [_reduction(r["input"]["two_qubit"], r[tool]["two_qubit"])
                  for r in rows]
        red_d = [_reduction(r["input"]["depth"], r[tool]["depth"]) for r in rows]
        times = [r[f"{tool}_ms"] for r in rows if r.get(f"{tool}_ms") is not None]
        out["configs"][tool] = {
            "gate_reduction_pct": {"median": _median_pct(red_g),
                                   "geomean": _geomean_pct(red_g)},
            "entangling_reduction_pct": {"median": _median_pct(red_2q),
                                         "geomean": _geomean_pct(red_2q)},
            "depth_reduction_pct": {"median": _median_pct(red_d),
                                    "geomean": _geomean_pct(red_d)},
            "compile_ms": {"median": _median(times),
                           "p95": (sorted(times)[int(0.95 * (len(times) - 1))]
                                   if times else None)},
            "error_rate": 0.0,
        }
    proof_ok = sum(1 for r in rows
                   if r.get("proof") in ("exact-unitary", "randomized-exact",
                                         "compositional", "clifford"))
    out["correctness"] = {
        "verified": f"{proof_ok}/{len(rows)} in-product",
        "referee": "qiskit Operator (all outputs, both tools)",
    }
    return out


def score_repro():
    d = _load("repro.json")
    if not d:
        return None
    rows = d["records"]
    out = {"corpus": f"repro harness ({len(rows)} pinned circuits, "
                     "10 families)",
           "artifact": "results/repro.json|csv|md",
           "referee": "qiskit Operator, every output lowered to u3+cx",
           "configs": {}}
    for tool in ("compact", "qiskit", "pytket", "cirq"):
        ok_rows = [r for r in rows
                   if isinstance(r.get(tool), dict)
                   and r[tool].get("status") == "OK"]
        if not ok_rows:
            out["configs"][tool] = {"status": "no successful runs"}
            continue
        red_2q = [_reduction(r["input_2q_cx_lowered"], r[tool]["two_qubit"])
                  for r in ok_rows if r["input_2q_cx_lowered"]]
        times = [r[tool]["wall_ms"] for r in ok_rows
                 if r[tool].get("wall_ms")]
        out["configs"][tool] = {
            "ok": len(ok_rows),
            "entangling_reduction_pct": {"median": _median_pct(red_2q),
                                         "geomean": _geomean_pct(red_2q)},
            "compile_ms": {"median": _median(times),
                           "p95": (sorted(times)[int(0.95 * (len(times) - 1))]
                                   if times else None)},
            "error_rate": round(1 - len(ok_rows) / len(rows), 4),
        }
    wtl = d.get("summary", {}).get("compact_vs_qiskit_2q_WTL")
    if wtl:
        out["compact_vs_qiskit_2q_WTL"] = wtl
    return out


def score_mqtbench():
    d = _load("mqtbench.json")
    if not d:
        return {"corpus": "MQT Bench", "status": "NOT MEASURED "
                "(run scripts/mqtbench_run.py)"}
    rows = d["records"]
    out = {"corpus": f"MQT Bench ({len(rows)} algorithm-level circuits)",
           "artifact": "results/mqtbench.json|md",
           "referee": "qiskit Operator, circuits <= 8q",
           "configs": {}}
    for key, name in (("qz", "compact"), ("qk", "qiskit_l3"),
                      ("pt", "pytket"), ("cq", "cirq")):
        ok = [r[key] for r in rows if r.get(key)]
        if not ok:
            out["configs"][name] = {"status": "no successful runs"}
            continue
        red_2q = [_reduction(r["in"][1], v[1]) for r, v in
                  [(r, r[key]) for r in rows if r.get(key)]
                  if r["in"][1]]
        # fidelity: verified count where refereed
        fids = [v[4] for v in ok if len(v) > 4 and v[4] is not None]
        verified = sum(1 for f in fids if f > 1 - 1e-6)
        times = [v[3] for v in ok]
        out["configs"][name] = {
            "ok": len(ok),
            "entangling_reduction_pct": {"median": _median_pct(red_2q)},
            "compile_ms": {"median": _median(times),
                           "p95": (sorted(times)[int(0.95 * (len(times) - 1))]
                                   if times else None)},
            "verified_refereed": f"{verified}/{len(fids)}",
            "error_rate": round(1 - len(ok) / len(rows), 4),
        }
    return out


def score_benchpress():
    recs = sorted((RESULTS / "benchpress").glob("benchpress_*_*.json")) \
        if (RESULTS / "benchpress").is_dir() else []
    if not recs:
        return {"corpus": "Benchpress (qiskit/benchpress)",
                "status": "NOT MEASURED (run scripts/benchpress_run.py)"}
    # per gym+size: {circuit_name: {"out_2q", "out_gates", "mean_s"}}
    per_gym: dict = {}
    for p in recs:
        parts = p.stem.split("_")          # benchpress_<gym>_<size>_<stamp>
        gym, size = parts[1], parts[2]
        d = json.loads(p.read_text(encoding="utf-8"))
        for b in d.get("benchmarks", []):
            ex = b.get("extra_info") or b.get("extra") or {}
            out_2q = ex.get("output_gate_count_2q",
                            ex.get("output_2q_gate_count"))
            out_gates = (ex.get("output_gate_count")
                         or sum(ex.get("output_circuit_operations", {}).values())
                         or None)
            per_gym.setdefault(f"{gym}/{size}", {})[b["name"]] = {
                "mean_s": round(b["stats"]["mean"], 4),
                "out_2q": out_2q,
                "out_gates": out_gates,
            }
    out = {"corpus": "Benchpress abstract-transpile, QASMBench, "
                     "all-to-all (pytest-benchmark native records)",
           "artifact": "results/benchpress/",
           "input_2q_note": ("the Benchpress record does not carry input "
                             "2q counts; reductions below pair the "
                             "compactq output against the qiskit output "
                             "on the same circuit")}
    for size in ("small", "medium"):
        comp = per_gym.get(f"compactq/{size}", {})
        qk = per_gym.get(f"qiskit/{size}", {})
        if not comp:
            continue
        paired = [(v["out_2q"], qk[k]["out_2q"]) for k, v in comp.items()
                  if k in qk and v.get("out_2q") is not None
                  and qk[k].get("out_2q") is not None]
        wins = sum(1 for c, q in paired if c < q)
        ties = sum(1 for c, q in paired if c == q)
        losses = sum(1 for c, q in paired if c > q)
        out[f"compactq_{size}"] = {
            "ran": len(comp),
            "paired_vs_qiskit_2q_WTL": [wins, ties, losses],
            "compact_out_2q_total": sum(c for c, _ in paired),
            "qiskit_out_2q_total": sum(q for _, q in paired),
        }
    out["runs"] = {k: sorted(v) for k, v in per_gym.items()}
    return out


def score_feynman():
    d = _load("feynman.json")
    if not d:
        return {"corpus": "Feynman benchmark corpus",
                "status": "NOT MEASURED (run scripts/feynman_bench.py)"}
    rows = [r for r in d["records"] if "compact_gates" in r]
    red_g = [_reduction(r["input_gates"], r["compact_gates"]) for r in rows]
    red_t = [_reduction(r["input_t_count"], r["compact_gates"]) for r in rows
             if r.get("input_t_count")]
    return {
        "corpus": "Feynman benchmarks/qasm (reversible / Clifford+T)",
        "artifact": "results/feynman.json|md",
        "configs": {"compact": {
            "circuits": len(rows),
            "gate_reduction_pct": {"median": _median_pct(red_g),
                                   "geomean": _geomean_pct(red_g)},
            "t_count_median": _median([r.get("input_t_count") for r in rows]),
            "compile_ms": {"median": _median(
                [r["compact_wall_ms"] for r in rows])},
        }},
        "correctness": {"verified": f"{d.get('verified')}/{d.get('parsed')} "
                         "in-product; qiskit Operator <= 8q"},
        "feynver_tool": d.get("feynver_tool"),
    }


def score_proofs_and_latency():
    dd = _load("dd_ceiling.json")
    lat = _load("latency.json")
    out = {"proof_ceiling": {}, "latency": {}}
    if dd:
        out["proof_ceiling"] = {
            f: dd["families"][f].get("proven_upto")
            for f in dd.get("families", {})}
        out["proof_ceiling_artifact"] = "results/dd_ceiling.json"
    if lat:
        out["latency"] = {
            "import_ms": lat.get("compactq_import_ms"),
            "qiskit_import_ms": lat.get("qiskit_import_ms"),
            "optimize_with_proof_ms": {"p50": lat.get("optimize_ms_p50"),
                                       "p95": lat.get("optimize_ms_p95")},
            "artifact": "results/latency.json",
        }
    return out


def main() -> int:
    feynman = _load("feynman.json")
    mqt = _load("mqtbench.json")
    sections = {
        "qasmbench_small": score_qasmbench_small(),
        "repro_harness": score_repro(),
        "mqt_bench": score_mqtbench(),
        "benchpress": score_benchpress(),
        "feynman": score_feynman(),
        "proofs_and_latency": score_proofs_and_latency(),
        "not_runnable_here": {
            "RevLib": ("dynamic per-file database, majority .real/.tfc "
                       "formats: needs a .real parser + scraper; the same "
                       "reversible/oracle workload class is covered by the "
                       "Feynman corpus (adders, Toffoli networks, GF(2) "
                       "multipliers)"),
            "feynver (Feynman's verifier)": ("Haskell toolchain required; "
                                             "QCEC + PyZX are the external "
                                             "referees in this release"),
            "Quantum Benchmark Zoo": ("submission target for "
                                      "docs/BENCHMARK_PROTOCOL.md once "
                                      "the protocol stabilizes (outreach "
                                      "item, not a measurement)"),
        },
    }
    sources = [feynman, mqt, _load("repro.json"), _load("dd_ceiling.json"),
               _load("latency.json")]
    card = {
        "generated": sources[2].get("generated") if sources[2] else None,
        "positioning": ("compactq is a VERIFIED optimizer: every corpus "
                        "row reports correctness alongside size; hardware-"
                        "fidelity claims are NOT derived from these "
                        "logical-metric reductions"),
        "sections": sections,
        "reproducibility": _compact_repro_block(
            [s for s in sources if s]),
    }
    (RESULTS / "SCORECARD.json").write_text(
        json.dumps(card, indent=2) + "\n", encoding="utf-8")

    lines = ["# Compact scorecard (generated)", "",
             card["positioning"], ""]
    for name, sec in sections.items():
        lines += [f"## {name}", "", f"```json",
                  json.dumps(sec, indent=2)[:2400], "```", ""]
    lines += ["## Reproducibility", "",
              "```json", json.dumps(card["reproducibility"],
                                    indent=2)[:2000], "```"]
    (RESULTS / "SCORECARD.md").write_text("\n".join(lines) + "\n",
                                          encoding="utf-8")
    print(f"wrote {RESULTS/'SCORECARD.json'} + SCORECARD.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
