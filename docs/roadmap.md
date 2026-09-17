# Roadmap

- **v0.1 (current)** — initial public release: the verified optimizer
  (peephole / commutation / CP / pure-Python KAK / Clifford-tableau /
  phase-polynomial passes, multi-pipeline search, exact + approximate
  modes, hardware `Target` objectives, SABRE-lite routing, exact
  placement, QASM2/3 IO, Qiskit bridge + plugin) and the open
  error-suppression stack (noise-aware layout, exact Pauli twirling,
  benefit-gated DD with a sequence family, MLE readout mitigation,
  exact ZNE, CDR, trajectory + stabilizer simulators, device metrics,
  classical shadows, MaxCut-QAOA solver, resource estimates, execution
  adapters)

## Next

- Clifford+T fragment resynthesis (CliffordSimp-class pass) for
  basis_trotter-class rings — the one measured 2q loss (240 vs 179).
  Investigated 2026-09-16 on the current codebase: the plateau at 240
  is confirmed across ALL current levers — optimize_search depth 3/4,
  winmerge merge_by_unitary (exact), kak_pass(force=True),
  search->optimize fixpoint, and the Clifford+T normal-form pass
  (inapplicable: the ring's majority angles are generic; only
  isolated sites are pi/4-quantized). Closing it needs a new
  fragment-resynthesis capability, not more tuning of existing passes.
- **Graph-form Clifford synthesis (ZX-class)** — the van den Nest
  decomposition U = R . CZ_Q . H_E . L: extract the stabilizer
  generators of U|0...0> from the tableau, GF(2)-reduce them to the
  graph Q (H-flips collected as the Hadamard layer), solve the basis
  bits from the signs, complete the residual local layer per wire
  against the 24-element single-qubit Clifford group, and prove the
  whole emission with clifford_equal.  Sparse graphs (rings,
  GraphState, QEC-prepare circuits) then re-synthesize with |E(Q)|
  CZ instead of CX chains.  Scoped 2026-09-16: the extraction core is
  specified; the residual local-layer bookkeeping (sign patterns per
  wire) is the remaining build.
- **MPS simulator (n ~ 100)** for 1D-local circuits: left-canonical
  MPS with SVD bond truncation, sequential Born sampling, numpy as a
  lazy optional import.  First scaffold rejected in review (sampling
  and environment bookkeeping must be rewritten carefully against
  stabsim/statevector as referees); scheduled as a dedicated build.
- **Pauli-propagation simulator (n ~ 100)** — Heisenberg stochastic
  propagation: CHP conjugation rules for the Clifford backbone,
  4-outcome unbiased branch sampling for non-Clifford rotations
  (cos/sin branch split), analytic depolarizing+readout damping on
  the observable.  Design finalized 2026-09-16; implementation
  deferred — the Y/H/S-gate sign rules must be implemented against
  the full CHP mod-4 phase algebra and refereed against the
  statevector at n <= 8 before n ~ 100 claims are made.
- **PEC (probabilistic error cancellation)** — attempted 2026-09-16:
  the quasi-probability inversion of the simulator's per-wire
  depolarizing channel is straightforward, but the sampled-operation
  bookkeeping for general Clifford gates requires the full
  Pauli-transfer-matrix inversion per gate; the scalar approximation
  fails unbiasedness at the minimal single-gate case (measured -0.95
  vs ideal -1).  Deferred until the PTM machinery is built (needs a
  per-gate PTM builder over the NoiseModel).
- Rust port of the remaining pass-engine hot loops (raw-latency crown)
- MQT Bench in CI (script exists: `scripts/mqtbench_run.py`)
- Hardware-backed suppression benchmarks (adapters shipped; needs
  credential-carrying runs)
- First PyPI publication (trusted-publisher release workflow ships;
  requires one-time PyPI + GitHub environment configuration)
- multi-OS native wheels via cibuildwheel (Windows wheels ship today;
  other platforms build from source with maturin)
