# Changelog

All notable changes to Q-PROOF Compact are documented here.

## [0.1.6] — 2026-09-18

### Fixed
- **Evidence integrity release.**  The 0.1.5 notes quoted an untracked
  6.0 s figure for the 256Q Ising/Trotter case while the committed
  gauntlet artifact recorded 146.8 s.  Both are real measurements with
  different instruments: the gauntlet times optimization **under the
  allocation tracker** (tracemalloc adds large overhead on
  allocation-heavy pipelines; 146.8 s), the profiling run without it
  completes in 6.0 s.  The README now reports both numbers with their
  measurement conditions.
- `compactq.verify` docstring said "tier 0-3 today" after tier 4
  shipped; corrected to 0-4 with the compositional-strategy caveat.
- `VERIFICATION_TIERS[4]` / docstring / tests fully consistent (T4
  shipped, T5 roadmap).

### Added
- `scripts/provenance.py` + embedded environment blocks in every
  benchmark artifact: git SHA, package version, dependency versions,
  Python, OS, CPU, timestamp.
- `scripts/check_provenance.py` — reproducibility gate validating
  provenance of every committed artifact (wired advisory into CI).
- Scalability report upgraded: runtime-tier classes (excellent /
  interactive / practical / batch / extreme) and P95 / max per width;
  the runtime ceiling is now the strict "practical" class (median
  ≤ 30 s), replacing a meaningless 25-hour threshold.
- `verify_segmented` exported at package top level
  (`compactq.verify_segmented`).

### Changed
- All benchmark artifacts regenerated on the 0.1.6 code at a single
  commit: bench_results, pytket referee, BQSKit head-to-head,
  scalability gauntlet.

## [0.1.5] — 2026-09-18

### Fixed
- **T4 API consistency**: `VERIFICATION_TIERS[4]` now reads
  `compositional_certificates` (shipped) — the module docstring, the
  tier table and the tests previously disagreed about whether T4 was
  roadmap or shipped.

### Added
- **T4.1 sequential-segment proof**: `compactq.compositional.
  verify_segmented(a, b, cuts_a, cuts_b)` — sound for ANY consecutive
  segmentations (per-segment phase factors commute with everything);
  each segment pair is dense-verified on its active qubits, so callers
  who can name narrow segments (Trotter layers, barrier boundaries,
  rewrite correspondence) get compositional proofs beyond the disjoint-
  block case.  Undecidable shapes return None, never a guess.
- **Heavy-circuit guardrail** in `optimize_search`: beyond 1,200 gates
  the expensive synthesis candidates (phase-polynomial re-synthesis,
  cross-pair KAK, unitary window merge, Cliffordize) are skipped and
  the cheap exact pipeline runs.  Measured effect: 256-qubit
  Ising/Trotter optimization drops from a **180 s gauntlet timeout to
  6.0 s** (-22% gates, never-grow holds), moving the measured runtime
  frontier through 256Q.  Circuits under the limit (all of QASMBench
  small) are behavior-identical (drift-tested).

### Changed
- **BQSKit head-to-head executed** (the review's largest evidence gap):
  BQSKit 1.2.1, its documented `bqskit.compile` pipeline, same
  normalized QASMBench small inputs, every output refereed —
  `results/bqskit.json`.  Result: **Compact 2q win 10 / tie 19 /
  lose 3**, geomean BQSKit/Compact 2q ratio 1.16; all 32 BQSKit outputs
  referee-clean.  Compact's measured losses (basis_change_n3, fredkin_n3,
  simon_n6) are documented.
- Scalability artifact regenerated on the 0.1.5 code with the
  heavy-circuit guardrail: **all three ceilings (correctness/stability,
  runtime, memory) now sit at 256Q** — the 256Q Ising/Trotter case
  dropped from a 180 s timeout to 6.0 s (median 34.3 s across its 8
  workloads), zero crashes/timeouts/incorrect outputs at any width.

## [0.1.4] — 2026-09-18

### Added
- **T4 compositional verification (prototype)**:
  `compactq.compositional` proves equivalence for block-structured
  circuits at ANY width — disjoint qubit components are dense-verified
  per block at each block's own width (soundness: block-diagonal
  unitaries factorize, so per-component proofs compose).  Integrated
  into `compactq.verify()` as tier 4; a 100-qubit circuit of 4-qubit
  blocks is now fully equivalence-proven without any 2^100 unitary.
  Tampered blocks are detected and the failing block's qubits are named.
- **Certificate format + `compactq verify` CLI**:
  `compactq.build_certificate()` emits the trust artifact — sha256
  hashes of both circuits' canonical QASM, before/after metrics
  (gates/2q/CX-equivalents/depth) and the verification evidence; and
  `compactq verify ORIGINAL.qasm OPTIMIZED.qasm -o certificate.json`
  re-checks the claim on another machine (exit 0 = equivalent).
- **Coupling-map presets** (`compactq.topology.coupling_preset`):
  line / grid / heavy-hex (IBM-style) / all-to-all for
  hardware-mapped routing benchmarks.
- Scalability gauntlet extended to **256 qubits**: the measured frontier
  is now reported as three different numbers — correctness/stability
  128Q (all workloads PASS), runtime boundary 128→256Q (256Q
  Ising/Trotter exceeds the 180 s per-case budget as a clean timeout,
  never a wrong answer), memory ≥ 256Q (peak 168 MB).  GHZ and
  random-Clifford remain **algebraically proven at 256 qubits**.
- README: "trust layer" positioning with the artifact architecture,
  verification tier table (T0-T5 with T4 shipped), certificate format
  docs, release ladder (0.2 Scalable+ / 0.3 Synthesis / 0.4 Hardware /
  0.5 Verified compiler / 1.0 Compact compiler), priority research
  order, fault-tolerant-compiler watch item.

## [0.1.3] — 2026-09-18

### Added
- **Scalability gauntlet** (`scripts/scale_gauntlet.py`): battle-tests
  widths 4 → 128 qubits × 8 workload families (GHZ, QFT, random
  Clifford, random general, Ising/Trotter, QPE-like,
  hardware-efficient ansatz, library random).  Every case runs in an
  isolated subprocess with a wall-clock timeout and tracemalloc peak
  memory; every output is verified at the strongest available tier —
  independent Qiskit Operator referee at small widths, exact algebraic
  tableau proof for Clifford workloads at ANY width, and the never-grow
  policy check (gates / 2q / depth / CX-equivalents) plus a
  determinism double-run everywhere.  Emits
  `results/scalability.json` + `results/scalability.md` with per-width
  PASS/TEST/FAIL and the three measured ceilings (correctness/stability,
  runtime, memory).
- README scalability section: the honest statement of measured ceilings
  (replacing guesses about "maximum supported width").

## [0.1.2] — 2026-09-17

### Added
- **`compactq.verify(original, optimized)`** — the public
  verification/evidence API: independently checks any two circuits and
  returns `{'equivalent', 'tier', 'method', 'global_phase_ignored',
  'runtime_ms'}`.  Tiers: 3 = algebraic Clifford-tableau proof (any
  width), 2 = dense unitary certificate (≤ 8q), 1 = randomized K-state
  sampling (numpy), 0 = prover unavailable; tiers 4 (compositional
  certificates) and 5 (formal proof) are reserved and documented.
- **Level-B metric**: `Circuit.cx_equivalent_count()` (SWAP = 3 CX);
  `bench_results.json` now emits CX-equivalents for Compact, Qiskit L3
  and the input baseline — the normalized 2-qubit cost alongside the
  logical count.
- **Committed referee artifacts**: `results/pytket_referee.json` — one
  auditable record per refereed pytket circuit (version, mode, output
  fidelity, status, metrics), making the inequivalence claims directly
  reproducible.
- `scripts/bqskit_bench.py` — the BQSKit column for benchmark suite v1
  (scaffold; produces `results/bqskit.json` once BQSKit is installed;
  no BQSKit numbers are claimed until then).
- README: 15-second evidence block, verification tier table, Level
  A/B/C metric definitions, "verification layer for quantum
  compilation" positioning, calibration-aware-compilation roadmap,
  ordered research sequence.

## [0.1.1] — 2026-09-17

### Added
- **Multi-objective optimization**: `objective=` on
  `optimize` / `optimize_deep` / `optimize_search` and `--objective` on
  the CLI.  Supported: `2q` (default), `depth`, `gate_count`, `latency`
  (depth alias), and `weighted` (minimize 1.0·2q + 0.1·depth +
  0.02·gates).  Each lexicographic objective reorders the same
  (2q, gates, depth) tuple — the chosen primary metric can never grow,
  by construction; the `2q` default is behavior-identical to 0.1.0
  (drift-tested).  Hardware-error weighting remains
  `compactq.target.optimize_for`.
- `results/` benchmark-artifact directory: `environment.json` (exact
  environment, commands and referee protocol of every documented run)
  and the artifact layout; `bench_results.json`/`bench_results.md`
  regenerated on the 0.1.1 code (commit-stamped, 43 rows).
- Objective-mode test groups (`tests/test_objectives.py`, 4 groups)
  wired into CI and the release pipeline.

### Changed
- Console script `compact` retired in favor of `compactq` on all
  platforms — the bare `compact` name collides with Windows' built-in
  `compact.exe` file-compression tool.  `compact-bench` and
  `python -m compactq` are unchanged.
- README: verification model stated as two exact tiers (algebraic /
  symbolic vs numerical unitary certificate, with the 1e-7 tolerance as
  part of the contract); "verified quantum compilation" positioning
  with an explicit fair-claims statement; capability matrix vs
  Qiskit/TKET/BQSKit/Cirq/Staq; suppression scope-of-claims note;
  roadmap elevated (Trotter/Hamiltonian simulation as the primary
  research target, scalable verification, benchmark suite v1).

### Removed
- Stale prototype-era files from the repository root: the v0.2 KAK
  engineering plan (`KAK_PLAN.md`, long since shipped) and an
  unreferenced pickle artifact (`best_before_clifford.pkl`).

## [0.1.0] — 2026-09-17

Initial public release. The public version series starts at 0.1.0; an
earlier private prototyping sprint (local iterations, 2026-09-10 to
2026-09-16) produced the engine and is retired — see the git history
for the full trail.

### Optimizer

- Peephole folding (1q runs → ≤3 canonical gates), commutation engine,
  SWAP templates, generalized diagonal sliding, cross-pair commutative
  window merging
- CP engine: `CX·RZ_t(θ)·CX = RZ_c(θ)·RZ_t(θ)·CP(−2θ)` — halves the CX
  count of QAOA/QFT/PEA-style phase ladders; CP merging; phase-polynomial
  re-synthesis of diagonal cores
- Pure-Python KAK/Weyl re-synthesis (Cross et al., arXiv:1811.12926):
  0/1/2/3-CX minimal templates by exact fidelity test, oracle-validated
  against Qiskit's Rust implementation (~1e-15 Weyl agreement)
- Clifford tableaux (Aaronson–Gottesman) with tableau proofs at any
  qubit count; Clifford+T normal-form pass; permutation-absorbing KAK
- Multi-pipeline `optimize_search()` with circuit-feature candidate
  selection; `optimize_large()` randomized K-state verification to
  ~30 qubits; approximate mode with measured per-block fidelity;
  hardware `Target` objectives (per-pair fidelities, native CZ/ECR/iSWAP
  rebase); SABRE-lite routing with error-weighted look-ahead; exact
  placement (provably minimal ≤6 wires); mcx/mcp parity networks
- OpenQASM 2.0 import/export, OpenQASM 3 export (+ common-subset
  import), Qiskit bridge + transpiler plugin, Cirq export
- Optional Rust kernels (`compactq-native`): KAK hot path, `sim_unitary`
  (raises the exact-proof ceiling to 8 qubits), `trace2`, parity-network
  BFS; auto-detected with silent pure-Python fallback

### Error suppression (open, local, exactly provable)

- `NoiseModel` (T1/T2, per-pair gate infidelity, readout confusion,
  quasi-static drift, durations; IBM calibration ingestion)
- Exact Pauli twirling (model-gated, portfolio mode), benefit-gated
  dynamical decoupling (xy4/xy8/xzx/pdd4 + auto), tensored readout
  mitigation with MLE (Richardson–Lucy EM), exact ZNE, CDR
- Zero-dependency trajectory simulator (n ≤ 14) and CHP stabilizer
  simulator (n ~ 60); device metrics (layer fidelity / EPLG);
  classical shadows; MaxCut-QAOA solver; resource estimates
  (T-count / T-depth); stim export; qiskit-runtime / Braket adapters
- `suppress_plan` / `suppress_execute` one-call pipeline with a
  per-stage `SuppressionReport` audit artifact; reachable from the
  CLI (`compact in.qasm --suppress`), the Python API, and the Qiskit
  plugin

### Project infrastructure

- 76-function zero-dependency test suite (`tests/run_tests.py`);
  1,018-check end-to-end gauntlet (`scripts/gauntlet.py`, 61 algorithm
  families × every public entry point) refereed by Qiskit's `Operator`
- Vendored QASMBench (pinned pnnl/QASMBench `357b942`); benchmark
  harnesses for QASMBench / MQT Bench / synthetic suites
- Version-consistency lint; CI and release workflows; MIT license;
  contributing / code-of-conduct / security policy

[0.1.0]: https://github.com/Q-PROOF/Compact/releases/tag/v0.1.0
