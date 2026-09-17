# Competitive research: Q-CTRL and Fire Opal

*Research date: 2026-09-15. Sources: q-ctrl.com, docs.q-ctrl.com, IBM/AWS quantum
docs, funding announcements. Compiled as the working brief for Q-PROOF Compact's
error-suppression module.*

## 1. The company: Q-CTRL

| | |
|---|---|
| HQ | Sydney, Australia (+ LA, Munich, Tokyo) |
| Founder/CEO | Prof. Michael J. Biercuk (quantum physicist, USyd) |
| Funding | ~US$133M+ total; Series B extended to US$113M (GP Bullhound, Oct 2024) — a quantum-software fundraising record |
| Products | **Fire Opal** (algorithm execution / error suppression), Boulder Opal (R&D-level quantum control), Black Opal (education), Fire Opal Vision |
| Positioning | "Quantum control" — software that makes unreliable qubits useful without changing hardware |
| Distribution | IBM Quantum (Qiskit Function "Performance Management"), IonQ cloud + Braket (Forte), Rigetti QCS (Ankaa-3), direct (fire.q-ctrl.com) |
| Model | Commercial SaaS, usage-based; closed source. No public pricing |

Headline demo records they market: 33q QPE, 75q verifiable entangled state,
103q constrained optimization (Network Rail), 127q QAOA on IBM hardware.

## 2. The target product: Fire Opal

**What it is:** a fully automated error-suppression pipeline sitting between the
user's algorithm and the QPU. One function call: user submits a circuit (or a
high-level solver problem), Fire Opal transpiles, suppresses, executes, and
returns corrected results. No knobs exposed in the standard flow.

**The technique stack (from their docs and publications):**

1. **Compiler / transpilation** — circuit-level optimization and basis
   translation; they claim validated outperformance of vendor compilers.
2. **Fidelity-aware layout** — qubit placement scored by calibration data
   (their "graph tooling" layout selection).
3. **Dynamical decoupling (DD)** — context-aware pulse sequences in idle
   windows to refocus dephasing and quasi-static drift.
4. **Pauli twirling / randomized compiling (RC)** — converts coherent gate
   errors into stochastic ones; averaged over randomized variants.
5. **Measurement error mitigation (MEM)** — readout-confusion calibration and
   inversion (tensored).
6. **Result safeguards** — failure detection, cost warnings, fewer shots.

They explicitly position **suppression > mitigation** in their marketing
(suppression acts during the run; mitigation post-processes).

**Claims:** 10x deeper circuits, 1000x accuracy and cost improvements, "random
to useful" outputs. Benchmarks are published as success-probability graphs
(Algorithmic benchmarking results page) comparing Fire Opal vs default Qiskit
execution on GHZ, BV, QAOA, QFT — always *their closed pipeline vs default*.

**The strategic weakness:** every claim is (a) closed-source and unauditable,
(b) tied to their paid cloud, (c) measured only against *default* execution —
nobody sees the technique stack, can certify exactness of the transformed
circuits, or run the pipeline locally on their own calibrations.

## 3. Our play: Q-PROOF Compact Suppression

Compact already owns the compiler layer (verified optimization, never-grow,
proof net). The suppression module completes a Fire-Opal-equivalent
technique stack as an **open, zero-dependency, exactly-provable library**:

| Fire Opal | Compact Suppression |
|---|---|
| Closed cloud service | Open-source library, PyPI, runs anywhere |
| Transformed circuits unauditable | Every pass exact up to global phase, proof-net verified |
| Requires their backend access | Works with any NoiseModel (dict ingestion, IBM calibration loading) |
| No knobs (or all knobs) | Sensible automation + full knobs |
| Paid, usage-priced | Free, local, data never leaves the machine |
| Suppression only on their stack | Composes with compactq optimizer AND user's own stack |

**Honesty line (non-negotiable):** Fire Opal's headline numbers come from real
hardware runs. Our benchmark evidence is simulator-based (density-matrix noise
simulation with coherent + stochastic + drift + decoherence + readout error
classes, seeds fixed, fully reproducible). The README says exactly that. We
compete on: verifiable exactness, zero dependencies, price, portability, and
demonstrated per-technique wins under controlled noise — not on out-claiming
their hardware results.

## 4. What "beat" means here (measurable)

1. **Every suppression layer must win in the error regime it targets**
   (twirl → coherent; DD → decoherence/drift; MEM → readout; layout →
   nonuniform calibration maps) — vs raw AND vs compiler-only.
2. **The full pipeline must beat raw execution in every scenario** (success
   probability factor > 1), and beat a qiskit-transpile-only baseline.
3. **Monotonicity by construction**: layers are applied only when the noise
   model predicts a benefit, so the pipeline can never *lose* fidelity —
   Fire Opal cannot offer this guarantee (no exactness proofs exist for it).
4. **One-call UX parity**: `suppress_execute()` matches their "single line of
   code" story without the cloud.

## 5. Repo decision

Build in this repo (`compactq` package), not a new one: suppression reuses the
circuit IR, proof net, hardware routing, noise model ingestion, and the PyPI
distribution (`compactq`) that already exist. A second repo would fork all of
that and dilute the brand. The suppression stack ships as `compactq.suppress*`
modules in the same package.

## 6. Sources

- [Fire Opal product page](https://q-ctrl.com/fire-opal)
- [Fire Opal docs — algorithmic benchmarking results](https://docs.q-ctrl.com/fire-opal/discover/adopt/algorithmic-benchmarking-results)
- [IBM Quantum — Q-CTRL Performance Management function](https://quantum.cloud.ibm.com/docs/en/guides/q-ctrl-performance-management)
- [AWS blog — Fire Opal error suppression for IonQ on Braket](https://aws.amazon.com/blogs/quantum-computing/improve-quantum-workload-performance-with-fire-opal-error-suppression-for-ionq-processors-on-amazon-braket/)
- [Q-CTRL Series B announcement](https://q-ctrl.com/blog/q-ctrl-sets-global-quantum-technology-fundraising-record-increasing-series-b-to-usd-113m-led-by-gp-bullhound)
- [IonQ × Q-CTRL partnership (Apr 2026)](https://thequantuminsider.com/2026/04/23/ionq-and-q-ctrl-partner-to-unlock-quantum-optimization-with-fire-opal-on-forte-quantum-processors/)
