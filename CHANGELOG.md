# Changelog

All notable changes to Q-PROOF Compact are documented here.

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
