# Error suppression

Compact's suppression stack applies the published techniques behind
commercial error-suppression services — randomized compiling, dynamical
decoupling, tensored measurement-error mitigation, fidelity-aware layout —
with compactq's distinguishing property: **every transformed circuit stays
exactly provable**, and every layer is gated so it can never be a net loss.

## One call

```python
from compactq import Circuit, Gate, suppress_execute, default_model

circ  = Circuit(4, [Gate("h", (), (0,))] +
                [Gate("cx", (), (j, j + 1)) for j in range(3)])
noise = default_model(4)
# noise = NoiseModel.from_qiskit_backend(backend)   # live IBM calibrations

result = suppress_execute(circ, noise)     # plan -> run -> mitigate
best   = max(result["probabilities"],
             key=result["probabilities"].get)
```

`suppress_execute` plans the pipeline, runs each variant through your
execution function (`run_fn(circuit, *, seed, shots) -> counts` — point it
at your hardware adapter) or the built-in zero-dependency simulator, pools
the counts, and applies tensored readout mitigation.

## The pipeline (and why the order is fixed)

1. **Optimize** — `optimize_search`, never-grow enforced.
2. **Expand** — CP entanglers become the exact CX form the passes handle
   (`expand_for_suppression`).
3. **Randomized compiling** — K Pauli-twirled variants per seed; adjacent
   same-wire Paulis are merged back into single gates (halves the pulse
   overhead). Each variant is unitarily identical to the input.
4. **Dynamical decoupling** — placed on *each variant's own schedule*
   (DD must see the final idle structure; Pauli merging must never compose
   DD sequences away). Benefit-gated per idle window: refocusable noise
   (quasi-static drift + dephasing, from the NoiseModel) must beat 3x the
   pulse overhead or nothing is inserted.
5. **Measurement mitigation** — tensored per-wire readout inversion on the
   pooled counts (quasi-probabilities, clipped and renormalized).

Every variant is proven equivalent to your input before use (dense proof
to 8 qubits, randomized-K beyond); anything doubtful degrades to the
original circuit.

## NoiseModel

```python
from compactq import NoiseModel

noise = NoiseModel.from_dict(5, {
    "t1_us": {0: 120, 1: 90, 2: 140},
    "t2_us": {0: 70, 1: 60, 2: 85},
    "readout": {0: [0.02, 0.03], 1: [0.01, 0.02]},
    "gate_infidelity": {"1q": 0.0004, "cx": 0.007},
    "drift_rate": {0: 0.0004, 1: 0.0006},   # rms rad/ns, what DD refocuses
    "durations_ns": {"1q": 40, "cx": 280},
})
```

`drift_rate` is the rms quasi-static Z rotation rate per wire — the
low-frequency noise dynamical decoupling exists to refocus; DD's benefit
gate reads the same field the simulator samples.

## Standalone layers

```python
from compactq import pauli_twirl, insert_dd, mitigate_counts, simulate_counts

twirled  = pauli_twirl(circ, seed=0)          # exact up to global phase
decoupled = insert_dd(circ, noise)            # benefit-gated XY4
counts   = simulate_counts(circ, noise, shots=1000, eps_coh=0.05)
probs    = mitigate_counts(counts, noise.readout, circ.num_qubits)
```

## Measured evidence

Two machine-readable, seed-fixed benchmarks (both gate-checked, exit
non-zero on regression):

- `python scripts/suppress_bench.py` → `suppress_results.json` — per-layer
  study over four error regimes; the full pipeline beats raw execution in
  **every** scenario (x1.00–2.72 by regime, largest on decoherence-dominated
  BV chains).
- `python scripts/head_to_head.py` → `head_to_head_results.json` — raw vs
  qiskit transpile(optimization_level=3) vs qiskit-O3+suppression vs
  compact-full; compact-full wins every scenario in aggregate and every
  cell vs raw.

Honesty line: these are simulator-based, fully reproducible measurements of
the published techniques — not hardware claims. Hardware runs plug the same
pipeline into your backend via `run_fn` and `NoiseModel.from_qiskit_backend`.

## Research positioning

See [fire_opal_research.md](fire_opal_research.md) for the competitive
analysis against Q-CTRL's Fire Opal that motivated this module.

## The completed pipeline

- **Device-space mapping**: `suppress_plan(circ, noise, coupling=...)`
  places on the calibration-watched subgraph (`layout_aware`) and routes
  with poison-sighted SABRE-lite (`route_aware`); every variant is proven
  equivalent to the device-space reference of the input, end to end.
- **SuppressionReport**: per-stage gates/2q/depth + proof level
  (exact-unitary / tableau / randomized-K / input-unchanged) + measured
  suppression factor; `result["report"]` prints as a table and
  round-trips through `to_dict()`.
- **MLE mitigation**: `mitigate_mle` (Richardson-Lucy EM) - always a
  physical distribution; `suppress_execute(..., mitigation="mle")`.
- **Execution adapters**: `qiskit_runtime(backend)` and
  `braket_device(device)` return `RunTarget(run_fn, noise_model)`.
- **Scale**: `compactq.stabsim` (CHP stabilizer simulation, n ~ 60) and
  `scripts/scale_bench.py` - aware+mem beats naive mapping at n=16..24
  (e.g. GHZ-16: 0.399 naive -> 0.500 aware -> 0.695 with mitigation).
- **ZNE**: `compactq.zne` - provable identity folding, Richardson /
  linear / quadratic / exponential extrapolation, one call via
  `zne_execute(circ, obs_wires, noise)`.

## Sequence families + full surface

- **DD sequences**: `insert_dd(..., sequence=)` and
  `suppress_plan(..., dd_sequence=)` accept `xy4`, `xy8`, `xzx`,
  `pdd4`, or `auto` (XY8 for windows >= 8 pulse durations).  Every
  sequence is a Pauli word with identity product - exact, verified
  per family.  Sequence gate: `python scripts/suppress_bench.py
  --dd-sequences` (every family beats raw everywhere; `auto` is
  best-or-tied in most cells).
- **CLI**: `compactq in.qasm --suppress --noise calibrations.json
  --report report.json -o suppressed.qasm --json` - the whole
  pipeline offline, with the per-stage proof report as an artifact.
- **Qiskit plugin**: `suppress_qiskit(qc, coupling_map=...,
  noise_model=backend)` and `make_suppression_pass(...)` for
  PassManager composition.
