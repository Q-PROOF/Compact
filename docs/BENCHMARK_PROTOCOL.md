# The Compact Verified-Optimizer Benchmark Protocol (v1)

An open, reusable benchmark specification for **verified** quantum
circuit optimizers.  Designed to be contributed to the Quantum
Benchmark Zoo once stable; immediately usable by anyone from a clean
clone of this repository.

## Principles

1. **Identical inputs.**  Every tool receives byte-identical circuits,
   generated or vendored with pinned digests (`results/MANIFEST.json`).
2. **One gate set before counting.**  All outputs are lowered to a
   common basis (`u3+cx`) before any metric is computed.
3. **Correctness is a metric, not an assumption.**  Every output is
   refereed by an independent implementation (Qiskit `Operator`,
   MQT QCEC, PyZX).  An inequivalent output invalidates the row.
4. **Multi-metric, never one number.**  Gate reduction (median AND
   geometric mean), entangling-gate reduction, depth reduction,
   compilation cost (median/p95/timeout rate), correctness
   (verified/total), per-config rows.
5. **Honest dispositions.**  A tool that cannot run is `NOT MEASURED`
   or `NOT RUNNABLE HERE` with the reason — never estimated.
6. **Reproducibility contract.**  Commit SHA, package versions, seeds,
   machine, and the exact command are embedded in every artifact.
7. **Proof reporting.**  For verified optimizers, the proof tier and
   (when applicable) certificate presence are part of the record.

## Corpora

| corpus | source | scope |
|---|---|---|
| QASMBench small/medium | vendored `third_party/QASMBench` (PNNL) | 43+ unitary cores |
| MQT Bench | `mqt-bench` pip package (TUM), algorithmic level | algorithmic circuits, ≤ ~30q |
| repro-harness suite | generated in-repo, digest-pinned | 10 families × 2 sizes |
| Benchpress | github.com/qiskit/benchpress, abstract-transpile workout | QASMBench via pytest-benchmark |
| Feynman | github.com/meamy/feynman `benchmarks/qasm` | 44 reversible / Clifford+T circuits |
| RevLib | revlib.org | planned: needs a `.real` parser + scraper |

## Result schema (per record)

```json
{
  "circuit": "<name>", "qubits": 6,
  "input":  {"gates": 100, "two_qubit": 40, "depth": 60, "t_count": 12},
  "config": "<tool + version + options>",
  "output": {"gates": 70, "two_qubit": 26, "depth": 44, "t_count": 8},
  "compile_ms": 12.3,
  "verified": true,
  "referee": {"tool": "qiskit.quantum_info.Operator", "version": "2.5.2",
               "fidelity": 1.0},
  "proof":  {"tier": 2, "method": "dd_full_unitary", "certificate": true}
}
```

## Runners (this repository)

| command | produces |
|---|---|
| `python scripts/run_benchmarks.py` | repro + DD-ceiling + latency |
| `python scripts/bench_json.py` | QASMBench small vs Qiskit L3 |
| `python scripts/mqtbench_run.py` | MQT Bench four-way |
| `python scripts/benchpress_run.py --gym compactq` (then `--gym qiskit`) | Benchpress native pytest-benchmark records |
| `python scripts/feynman_bench.py` | Feynman corpus verified optimization |
| `python scripts/scorecard.py` | `results/SCORECARD.json|.md` aggregate |

## Claim policy

Corpus-specific statements only — e.g. "on the pinned QASMBench-small
corpus at the shared basis, Compact reduced median 2q count by X%,
refereed-equivalent on Y/Z circuits, median wall T ms".  Percentage
reduction is `100·(in−out)/in`.  Logical-metric reductions are compiler
metrics and are never presented as hardware-fidelity improvements.
Losses are reported with the same prominence as wins.
