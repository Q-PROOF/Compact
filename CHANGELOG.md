# Changelog

All notable changes to Q-PROOF Compact are documented here.

## [Unreleased]

### Changed
- Console script `compact` retired in favor of `compactq` on all
  platforms — the bare `compact` name collides with Windows' built-in
  `compact.exe` file-compression tool. `compact-bench` and
  `python -m compactq` are unchanged.

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
