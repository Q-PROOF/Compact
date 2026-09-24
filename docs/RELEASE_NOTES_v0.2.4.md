# Q-PROOF Compact v0.2.4 — reproducible verified circuit-optimization results

**Tagline:** smaller circuits, proven — and now every number in this
repository regenerates from one command against a pinned manifest, on
standard corpora (QASMBench, MQT Bench, Benchpress, the Feynman
benchmark corpus), with the proof tier of every output attached.

## Highlights

- **Reproducibility contract**: the benchmark suite is pinned by
  `results/MANIFEST.json` (per-circuit QASM digests + referee/tool
  versions); `python scripts/run_benchmarks.py` regenerates the
  artifact families and aborts loudly on environment drift.
- **Four benchmark corpora, results committed**: QASMBench small
  (43 circuits), MQT Bench (algorithm-level, four-way),
  [Benchpress](https://github.com/qiskit/benchpress) abstract-transpile
  workout (native pytest-benchmark records via a new `compactq` gym in
  `benchpress_integration/`), and the Feynman benchmark corpus (44
  reversible/Clifford+T circuits — adders, Barenco Toffoli networks,
  GF(2) multipliers).  Aggregate multi-metric scorecard:
  `results/SCORECARD.md`.
- **Correctness fixes found by the pre-release red-team** (all pinned
  by `tests/test_adversarial_v024.py`):
  - the DD prover returned a false INEQUIVALENT verdict on ≥7-control
    MCX self-comparisons — fixed by max-weight (QMDD-standard)
    normalization + a self-consistency norm guard; QFT-family DD
    builds ~7× faster; Grover-MCX exact ceiling 10q → 12q.
  - `cx q[0], q[0]` (repeated-wire gates) were silently accepted by
    the IR — now rejected at construction.
  - qelib1 `ch` failed QASM import — now imported exactly.
  - `verify()` had no time/memory bounds past the DD tier — the
    prover cascade now declines on wall-clock deadlines
    (`COMPACTQ_DD_DEADLINE_S`, `COMPACTQ_VERIFY_DEADLINE_S`) and a
    tunable node budget (`COMPACTQ_DD_MAX_NODES`) instead of grinding
    or thrashing; declined verdicts are loud (`None` / tier 0), never
    guesses.
- **Use Compact inside your toolchain**: the Qiskit
  `TransformationPass` (fixed — the previous facade crashed
  PassManager) and a new pytket pass, both preserving the verification
  verdict as inspectable pass output; loss-free pytket bridge.
  `pip install compactq[qiskit]` / `compactq[pytket]`.
- **Governance docs**: VERIFICATION.md (tier-by-tier proof
  methodology, every boundary backed by a named test),
  CORRECTNESS.md (adversarial-testing invitation + found-bugs ledger),
  OUTREACH.md (independent-validation invitations — tracked open),
  compactq-check/INDEPENDENCE.md, docs/BENCHMARK_PROTOCOL.md (the open
  protocol intended for Quantum Benchmark Zoo submission).

## Honest scope

Not claimed this release: independent validation (invitations sent /
tracked in OUTREACH.md), RevLib as an automated corpus (needs a
`.real` parser; its reversible workload class is covered by the
Feynman corpus), `feynver` as an external verifier (Haskell
toolchain; QCEC + PyZX serve as external referees).  QFT-family
exact-DD proofs remain at 18q on the default budget — the Rust DD
kernel is still the roadmap item for 20-25q there.

## Numbers (2026-09-23 runs, commit-pinned artifacts)

See `results/SCORECARD.md` for the full multi-metric table.
Headline rows (identical inputs, u3+cx-normalized counting, Qiskit
`Operator` refereed):

- repro harness (20 pinned circuits, 10 families): compact mean 2q cut
  10.2%, median wall 22 ms (proof included) vs Qiskit L3 14.6% / 7 ms;
  head-to-head 2 W / 12 T / 1 L.
- Feynman corpus: verified optimization with large reductions on
  reversible workloads (e.g. GF(2) multipliers: gf2^128_mult 279,419 →
  213,881 gates; barenco_tof_5 218 → 122).
- DD proof ceiling: QFT 18q exact (declines at 20q), Grover-MCX 12q,
  Trotter/GHZ through 30q — with before/after deltas in
  `results/dd_ceiling.json`.

**Install:** `pip install compactq` (zero dependencies, Python ≥ 3.9).
Extras: `compactq[qiskit]`, `compactq[pytket]`, `compactq[bench]`.
