<p align="center">
  <img src="docs/assets/banner.svg" alt="Q-PROOF Compact — the verified quantum circuit optimizer. Smaller circuits, proven." width="100%">
</p>

<p align="center">
  <a href="https://pypi.org/project/compactq/"><img alt="PyPI version" src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fpypi.org%2Fpypi%2Fcompactq%2Fjson&query=%24.info.version&label=pypi&color=orange"></a>
  <a href="https://pypi.org/project/compactq/"><img alt="Python" src="https://img.shields.io/badge/python-%E2%89%A5%203.9-blue"></a>
  <a href="https://pypistats.org/packages/compactq"><img alt="Downloads" src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fpypistats.org%2Fapi%2Fpackages%2Fcompactq%2Frecent&query=%24.data.last_month&label=downloads&suffix=%2F+month&color=blue"></a>
  <a href="https://github.com/Q-PROOF/Compact/releases"><img alt="Release" src="https://img.shields.io/github/v/release/Q-PROOF/Compact"></a>
  <a href="https://github.com/Q-PROOF/Compact/blob/main/LICENSE"><img alt="License" src="https://img.shields.io/pypi/l/compactq"></a>
  <a href="https://github.com/Q-PROOF/Compact/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/Q-PROOF/Compact/actions/workflows/ci.yml/badge.svg"></a>
</p>

**The verified quantum circuit optimizer.** Smaller circuits, proven. Compact takes a quantum circuit and returns an equivalent one that is smaller and shallower — with verification attached to every answer: exact algebraic proofs where the circuit structure allows them, and a numerical whole-unitary certificate otherwise. On any doubt, your input is returned unchanged.

**Current evidence** — 43 QASMBench circuits + 44 MQT Bench circuits refereed
by Qiskit's `Operator`; a 1,018-check end-to-end gauntlet; 84 test functions;
**0 incorrect Compact outputs** in any refereed set; measured Qiskit L3 /
pytket / Cirq comparisons; hardware-error objectives, SABRE-lite routing,
and an open error-suppression stack.
**Current limitation:** BQSKit and hardware-mapped benchmarking are still
being expanded (`results/`, roadmap).

```python
import compactq
from compactq.benchmarks import qft

opt = compactq.optimize_search(qft(4))
print(opt.stats())
```

`pip install compactq` — no dependencies, no account, no cloud. Pure Python ≥ 3.9.
(From source: `pip install git+https://github.com/Q-PROOF/Compact.git`.)

## Why

Every gate you remove from a quantum circuit removes noise. Compact takes a circuit and
returns an **equivalent** one that is smaller (fewer gates), shallower (lower depth),
with priority on cutting 2-qubit gates — the dominant error source on today's hardware.
The objective is exactly lexicographic: **(2-qubit count, total gates, depth)** —
compactq accepts a slightly deeper circuit when it removes gates, never a larger
2-qubit count.

## Correctness model (the part that matters)

- Every rewrite is **exact**: the unitary is preserved up to global phase.
- For circuits ≤ 8 qubits (6 without the optional native kernel), `optimize_search()`
  **proves** equivalence via full-unitary comparison before returning; on any numerical
  doubt it returns your circuit unchanged. (`verify=False` skips the proof for large
  circuits — and is itself fuzzed against a Qiskit referee up to 10 qubits.)
- The KAK/Weyl re-synthesiser is cross-validated against Qiskit's Rust implementation
  (Weyl coordinates agree to ~1e-15 on randomized SU(4)s) and every re-synthesized
  block is re-verified against its own 4×4 unitary before it can replace anything.
- Clifford blocks are re-synthesized with **tableau proofs**, exact at any qubit count.
- The test suite (76 test functions in `python tests/run_tests.py`, plus
  objective-mode tests in `tests/test_objectives.py`) includes property tests over thousands of random circuits
  and a regression suite for the classic optimizer bugs (gate-order reversals,
  reversed-CX "cancellations", Euler-angle wrapping, ZZ-identity sign errors).
  `scripts/gauntlet.py` adds a 1,018-check end-to-end gauntlet: 61 realistic algorithm
  families × every public entry point, QASMBench through compactq's own importer, >8q
  no-verify soundness and CLI/bridge/round-trip checks — all refereed by Qiskit's
  `Operator`, never by compactq's own math.

**What "verified" means — two exact tiers, stated precisely:**

- **Algebraic / symbolic (exact, no floating point).** Clifford blocks are
  proven with stabilizer tableaux (`clifford_equal`) at any qubit count;
  every rewrite (peephole foldings, the `CX·RZ·CX → RZ·RZ·CP` identity,
  SWAP templates, Euler refolds) is an algebraic identity checked at each
  application; every KAK re-synthesis is re-verified against its own 4×4
  block unitary before it can replace anything.
- **Numerical unitary comparison (general circuits ≤ 8 qubits).** The
  optimized circuit is accepted only when the Hilbert–Schmidt fidelity
  |Tr(U†V)|/d against the input exceeds 1 − 1e-7.  This is a machine-checked
  numerical certificate — strong, but not a formal proof-assistant-grade
  certificate; the tolerance is part of the contract, and on any numerical
  doubt the original circuit is returned unchanged.

The public API exposes this as grades — `compactq.verify(a, b)` returns the
verdict **and** the evidence tier:

```python
>>> compactq.verify(ghz(5), ghz(5))
{'equivalent': True, 'tier': 3, 'method': 'clifford_tableau',
 'global_phase_ignored': True, 'runtime_ms': 1}
```

| tier | evidence grade | status |
|---|---|---|
| 0 | none / prover unavailable | today |
| 1 | randomized state sampling | today (numpy) |
| 2 | full-unitary numerical equivalence | today (≤ 8q) |
| 3 | local algebraic certificate (tableau) | today (Clifford) |
| 4 | compositional certificates | roadmap |
| 5 | formal proof | roadmap |

## What's inside

- **Peephole folding** — maximal 1-qubit runs resynthesized to ≤3 canonical gates
  (named-Clifford recognition, single-RX/RY recovery, exact H·P / P·H two-gate forms,
  RZ-RY-RZ Euler), with a final single-`u3` fold for gate-count polish.
- **Commutation engine** — self-inverse CX/CZ cancellation across provably-commuting
  gates, X-on-target / diagonal-on-control slides, SWAP templates, generalized
  diagonal sliding, cross-pair commutative window merging.
- **CP engine** — `CX·RZ_t(θ)·CX = RZ_c(θ)·RZ_t(θ)·CP(−2θ)`: exact ZZ-phase
  extraction that halves the CX count of QAOA/QFT/PEA-style phase ladders; adjacent
  CP merging; phase-polynomial re-synthesis of diagonal cores.
- **Pure-Python KAK/Weyl synthesis** (`compactq/kak.py`) — the numerically-stable
  simultaneous-diagonalization algorithm (Cross et al., arXiv:1811.12926 App. B):
  real-symmetric Jacobi eigensolver of Re/Im magic-basis parts, Weyl-chamber
  canonicalization, and minimal-CX circuit templates (0/1/2/3 CX by exact
  fidelity test). No numpy, no qiskit, no pytket needed in the core.
- **Clifford stabilizer tableaux** (`compactq/stabilizer.py`) — Aaronson-Gottesman
  CHP with phase-exact Pauli conjugation (single-qubit gates conjugated
  numerically, so a wrong sign convention cannot ship silently) and AG block
  resynthesis with GF(2) sign correction. `is_clifford` and `clifford_equal`
  give exact Clifford equality at any qubit count without building 2ⁿ unitaries.
- **Approximate mode** (`compactq.approximate`, `compactq.target`) — 2q blocks re-synthesized
  into cheaper CX classes with *measured* per-block fidelity guarantees; hardware-aware
  `Target` objective (per-pair CX fidelities, CX-direction flipping) and greedy
  fidelity-budget allocation. Reported fidelities include the approximation cost.
- **Hardware-aware routing** (`compactq/hardware.py`) — SABRE-lite SWAP insertion with
  error-weighted look-ahead, optional permutation restore, `rz-sx-x` 1q translation.
- **Multi-controlled gates** — `mcx`/`mcp` expanded via parity networks; the Qiskit
  bridge boundary-decomposes anything else a `QuantumCircuit` may carry (`mcphase`,
  `ccx`, `ecr`, `iswap`, `cu`, `rxx`, ...) loss-free into the supported basis.
- **Verification net** — full-unitary fidelity proofs for ≤8-qubit circuits.
- **Optional Rust kernels** (`native/`) — a PyO3 wheel (`compactq-native`, abi3
  stable ABI, Python ≥3.9) that accelerates the KAK hot path (block unitaries,
  determinants) and raises the exact-proof ceiling to 8 qubits (`sim_unitary`).
  Windows wheels ship with each release; other platforms build from source with
  `maturin` (see Development). compactq auto-detects the kernel and silently
  falls back to the pure-Python path — the zero-dependency contract never changes.

## Positioning: verified quantum compilation

Compact's thesis is not "another optimizer" — it is **verified quantum
compilation**: don't just optimize the circuit, return evidence that the
transformation preserved the computation.

```
              quantum program
                    │
                    ▼
            ┌───────────────┐
            │   Compact     │
            │  optimization │
            └───────┬───────┘
                    │
        ┌───────────┼───────────┐
        ▼           ▼           ▼
    synthesis    routing    suppression
        │           │           │
        └───────────┼───────────┘
                    ▼
        verification / certificate
                    │
                    ▼
            optimized circuit
```

**Fair claims.** Measured, independently refereed evidence supports
competitive logical optimization — 2-qubit reduction, depth, structured
families — versus Qiskit L3, pytket and Cirq on the suites below.  Not yet
established, and not claimed: blanket superiority over the full
Qiskit/TKET/BQSKit compiler stack (the Benchpress study, *Nature
Computational Science* 2025, found no single tool dominating across
construction, compilation and topology), advantage on hardware-mapped
circuits, real-device suppression advantage, and formal
proof-assistant-grade certificates.  Closing those gaps is the roadmap.

**The strategic frame:** Compact is the **verification layer for quantum
compilation**.  Qiskit, TKET and BQSKit are superb compilers; Compact is
complementary — it optimizes, and it independently verifies anyone's
output:

```python
qiskit_out = transpile(circ, optimization_level=3)
compactq.verify(circ, compactq.from_qiskit(qiskit_out))
# {'equivalent': True, 'tier': 2, 'method': 'full_unitary', ...}
```

Compile with anything; trust the circuit only when the evidence tier says so.

## How Compact competes

Three axes, all measured on identical inputs (QASMBench unitary cores, level-0
normalized, same basis, 2026-09-08 run unless noted):

Capability matrix (✅ shipped · ◐ partial · ❌ not claimed):

| capability | Compact | Qiskit | TKET | BQSKit | Cirq | Staq |
|---|---|---|---|---|---|---|
| logical optimization | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| 2-qubit unitary synthesis | ✅ | ✅ | ✅ | **strong** | ✅ | ✅ |
| layout / routing | ✅ | ✅ | ✅ | ✅ | ✅ | ◐ |
| hardware-aware objective | ✅ | ✅ | ✅ | ✅ | ◐ | ◐ |
| approximate optimization | ✅ | ◐ | ◐ | ✅ | ◐ | varies |
| in-product output verification | **✅** | external | external | external | external | external |
| measured head-to-head vs Compact | — | ✅ below | ✅ below | not yet | ✅ (MQT) | not measurable here |

(The cross-compiler reference point is the Benchpress suite — Qiskit, TKET,
BQSKit, Cirq, Staq and others; adding a measured BQSKit column is a named
roadmap item.)

**Metric levels.** A raw "2-qubit gate count" is a logical metric, not a
hardware-cost metric — a SWAP is three CX, and CP/RZZ/ECR costs are
target-dependent.  Benchmarks therefore report three levels:

- **Level A — logical:** 2Q operations, total operations, depth
  (`two_qubit_count`, `gate_counts`, `depth`).
- **Level B — normalized:** CX-equivalent 2Q count (SWAP = 3 CX) —
  `Circuit.cx_equivalent_count()`, emitted for Compact **and** Qiskit in
  `bench_results.json` (Level B for the remaining competitors lands with
  benchmark suite v1).
- **Level C — hardware:** native 2Q gates, weighted error cost, expected
  fidelity, estimated duration — via `compactq.target`, and the planned
  MQT four-level (algorithmic / target-independent / native / mapped)
  benchmark suite v1.

**1. Optimization quality** — win-or-tie on 2-qubit count vs Qiskit L3 on 30/32
QASMBench small circuits (11 wins / 19 ties / 2 losses). One loss
(`basis_test_n4`, 12 vs "6") is an accounting artifact: Qiskit L3 *elides* 2 SWAP
gates into its layout metadata instead of the circuit — reified, they cost 6
CX-equivalents, i.e. parity with compactq's 12; the other (`basis_trotter_n4`,
240 vs 179) is a genuine open case, on the roadmap. Vs pytket
`FullPeepholeOptimise`, see axis 2: its headline 2q
numbers on small kernels come from invalid outputs. On freshly generated
**MQT Bench** algorithm circuits (44 circuits, every output refereed, table
below) compactq **wins 25 / ties 19 / loses 0** vs Qiskit L3 — QAOA 2q cut in
half at every size (8→4, 28→14, 56→28), QFT 14→8, QPE 9→5, Grover
52→44 — and **25 / 19 / 0** vs Cirq's `CZTargetGateset` optimizer. The
v1.3 CX-phase rewrite (matched CX pairs around folded T/phase runs collapse
to CP) closes every previously-lost 2q account: adder_n10 65→**57**,
adder_n4 10→**7**, toffoli 6→**5**, cdkm adders **15/29** vs pytket's
16/31, Grover −16%.

**2. Precision** — Compact *proves* every output before returning (full-unitary
fidelity proof ≤ 8q, tableau proofs for Clifford blocks at any size; it returns
your input unchanged rather than ship an unprovable result).
An independent referee (Qiskit `Operator`, |Tr(A†B)|/d) over 29 QASMBench circuits:

| tool | inequivalent outputs | what the user must do to be safe |
|---|---|---|
| **Compact** | **0 / 29** (proof attached in-product) | nothing |
| Qiskit 2.5 L3 | 1 / 29 *as shipped*: `basis_trotter_n4` returns fid 0.25 because 2 SWAPs sit in `qc.layout` metadata, not the circuit | inspect `qc.layout` and re-apply elided permutations |
| pytket `FullPeepholeOptimise()` (default) | **6 / 30 wrong circuits** (fresh 2.18.1 referee: basis_test_n4 fid 0.50, grover_n2 0.50, hs4_n4 0.25, iswap_n2 0.50, qec_en_n5 0.25, sat_n7 0.016); `replace_implicit_wire_swaps()` does **not** repair them | run in `allow_swaps=False` mode — but 2.18.1's safe mode is itself inequivalent on 2/30 (basis_trotter_n4, sat_n7), and its valid outputs still never beat compactq (9 / 19 / 0) |
| pytket safe mode | 2 / 30 (basis_trotter_n4 fid 0.5, sat_n7 fid 0.016) | re-referee before trusting |

**3. Latency** — the honest metric is the *same job*: an optimized circuit you can
trust. compactq's price includes the proof; competitors need an external
`Operator`-fidelity check afterwards (which is exactly what caught the failures
above).

| same job (optimize + verify) | compactq (proof included) | Qiskit L3 raw | Qiskit + verify | compactq wins | pytket + verify (compactq wins) |
|---|---|---|---|---|---|
| QASMBench, 29 circuits | **1.86 s** | 0.52 s | 0.88 s | **18/29** | 11.0 s, **29/29** |
| MQT Bench, 44 circuits | **4.88 s** | 0.85 s | 1.97 s | 22/44 | — |

(timings from the 2026-09-08 instrumentation run; the correctness side was
re-verified fresh on 2026-09-16/17 — the full gauntlet passes 1018/1018 on the
current code, and the referee counts above are from that run.)

compactq wins the same-job race on the **majority of circuits** (18/29 QASMBench,
22/44 MQT — tally per the table above) and is always 4–6x faster than pytket;
Qiskit's Rust pass engine keeps the raw-speed crown, especially on large circuits
(closing that gap needs compactq's pass engine itself in Rust — on the roadmap).
End-to-end optimizer+prover throughput doubled across the prototyping phase: the
1,018-check gauntlet runs 150 s → 75 s and the trotter6 search 4.9x (1.18 s →
0.24 s), via native-kernel verification at all qubit counts, a value-keyed
unitary memo, fused 1q-matrix comparison, memoized gate matrices, and the Rust
`trace2` fidelity kernel.

Competitor landscape (measured columns from the runs above; others qualitative):

| tool | scope | output verification | 2q optimization | notes |
|---|---|---|---|---|
| **Compact 0.1 (Q-PROOF)** | logical optimization, exact + approximate | **yes — in-product proof** | **best measured** (30/32 win-or-tie QASMBench; 25/19/0 vs Qiskit and Cirq on MQT; 9/19/0 vs pytket-2.18 safe mode) | same-job latency won on the majority of circuits; zero-dependency core, optional Rust |
| Qiskit 2.5 transpiler (L3) | full transpilation stack (layout/routing/noise-adaptive) | none | strong; parity with compactq once permutations are reified | Rust-fast, huge ecosystem |
| pytket 2.18 (`FullPeepholeOptimise`) | logical optimization | none — **6/30 wrong by default** (fresh referee; safe mode now also 2/30 wrong: basis_trotter_n4, sat_n7) | on its 24 valid default-mode outputs: compactq wins 8 / ties 15 / loses 1 (basis_trotter_n4); safe mode: compactq **9 / 19 / 0** | the SWAP-elision pitfall persists in 2.18.1 |
| Cirq 1.7 (`optimize_for_target_gateset`) | construction + target-gateset optimization | none — but its optimizer is exact on our whole MQT run | weak: 23 losses / 0 wins vs compactq on MQT; grows some circuits | fast (tens of ms) |
| staq (softwareqinc) | synthesis/optimization toolchain | n/a | n/a | **not measurable in this environment**: no PyPI distribution (the PyPI `staq` package is an unrelated C decompiler) and building it needs a C++17 toolchain that is absent here |
| MQT Bench (mqt-bench, ALG level) | benchmark *generator* (44 circuits run here) | — | — | integrated as a first-class suite: `python scripts/mqtbench_run.py` |

**DD sequences, every surface** — suppression goes deeper and ships everywhere:
a dynamical-decoupling **sequence family** (`DD_SEQUENCES`: xy4, xy8,
xzx, pdd4) with `auto` selection by window length (`dd_sequence=`
through `suppress_plan`/`suppress_execute`, sequence gate
`--dd-sequences`); the pipeline is reachable from **every surface** —
`compactq in.qasm --suppress [--noise noise.json] [--report
report.json] [--json]` on the CLI, and `make_suppression_pass()` /
`suppress_qiskit()` in the qiskit plugin (calibration ingestion from a
BackendV2, coupling-map aware, proof-verified variant 0).

**End-to-end pipeline** — noise-aware layout +
SABRE-lite routing wired into `suppress_plan` (device-space mapping with
per-edge calibration sight), a `SuppressionReport` artifact (per-stage
gate/depth/proof-level records + measured suppression factor), MLE
measurement mitigation (Richardson-Lucy EM - always a physical
distribution, beats clipped inversion 3-5x on injected confusion), thin
execution adapters (`compactq.adapters.qiskit_runtime` /
`braket_device`), model-gated twirling (`coherent_fraction`) with
portfolio `twirl_fraction`, an optimization regression gate
(`scripts/bench_gate.py`), a CHP stabilizer scale simulator
(`compactq.stabsim`, n~60; Y-convention audit discharged by the fuzz
suite) with a gate-checked n=16..24 scale benchmark
(`scripts/scale_bench.py`), and exact zero-noise extrapolation
(`compactq.zne` - provable identity folding + Richardson/poly fits).

## Error suppression — the open stack, one call

Fire Opal (Q-CTRL) sells automated error suppression as a closed cloud
service: transpilation, fidelity-aware layout, dynamical decoupling, Pauli
twirling, measurement mitigation — one function call, no knobs, results
that cannot be audited. Compact ships the same technique stack as an open,
zero-dependency, **exactly provable** library that runs locally against
any noise model:

| Fire Opal | Compact suppression |
|---|---|
| closed cloud, paid usage | open source, free, local — data never leaves the machine |
| transformed circuits unauditable | every pass exact up to global phase, proof-net verified before use |
| only on their supported backends | any `NoiseModel` (plain dict, or live IBM calibration loading) |
| black-box pipeline | automated default + every layer available standalone |
| suppression only on their stack | composes with any compiler — works on qiskit-O3 output too |

```python
from compactq import Circuit, Gate, suppress_execute, default_model

circ  = Circuit(4, [Gate("h", (), (0,))] + [Gate("cx", (), (j, j+1)) for j in range(3)])
noise = default_model(4)                      # or NoiseModel.from_qiskit_backend(backend)
result = suppress_execute(circ, noise)        # plan -> run -> mitigate
best   = max(result["probabilities"], key=result["probabilities"].get)
```

Pipeline (fixed order, each layer exact and gated): optimize → expand
untwirlable CP entanglers → K Pauli-twirled variants → dynamical decoupling
placed on each variant's own schedule → tensored readout mitigation.
DD is **benefit-gated**: a window is decoupled only when the model's
refocusable noise (quasi-static drift + dephasing) beats 3x the pulse
overhead, so the pass can never be a net loss. Every variant is proven
equivalent to your input before it is used.

**Measured** (density-matrix noise simulation — coherent overrotation per
2q gate with per-site axes, gate depolarizing, T1/T2 with ASAP scheduling,
quasi-static drift, readout confusion; fixed seeds; reproducible with
`python scripts/suppress_bench.py` and `python scripts/head_to_head.py`):

- full pipeline vs raw execution: wins in **every** scenario — success
  probability x1.00–1.01 (coherent-dominated), x1.04–**2.72** (decoherence),
  x1.17–1.31 (readout-dominated), x1.06–1.13 (combined)
- head-to-head (`head_to_head_results.json`): raw vs qiskit-O3 vs
  qiskit-O3+suppression vs compact-full — mean success probability
  coherent **0.825** (qiskit-O3 0.822, raw 0.808), decoherence **0.660**
  (0.431, 0.418), combined **0.766** (0.688, 0.677); single best cell
  BV n=6 decoherence **0.729 vs 0.138** (5.3x). The suppression stack
  also lifts qiskit-O3 output when composed onto it.
- one-cell honesty note: on coherent-dominated QAOA n=4, bare qiskit-O3
  beats *every* suppressed pipeline (including qiskit-O3+suppression) —
  randomized compiling trades coherent-error cancellation for stochastic
  robustness; the aggregate still favors the pipeline.

The zero-dependency trajectory simulator (`compactq.simulate_counts`)
reproduces the same physics for n <= 14 without numpy, so the whole
pipeline runs — and is tested — anywhere Python runs.

**Honesty line:** Fire Opal's headline numbers come from real hardware;
ours are simulator-based by construction (the techniques are the published
ones — Mundada et al. randomized compiling, XY4 decoupling, tensored
readout inversion). What we add is what no closed service can offer:
per-layer exactness proofs, benefit gating with a no-net-loss argument,
and a protocol anyone can re-run.

**Scope of claims.** Suppression competes in a different ecosystem (Q-CTRL
Fire Opal, Mitiq, Qiskit Runtime, Superstaq), where the correct benchmark is
not gate count but: measured fidelity improvement over raw execution,
sampling overhead, runtime overhead, and estimator bias/variance — on
experimental noise models and, ultimately, hardware.  Compact ships the
technique stack with simulator-based evidence only; hardware-backed
suppression numbers are explicitly *not* claimed (execution adapters are
shipped; credential-carrying runs are pending).

All numbers below are from a single re-run (2026-09-16, Python 3.11.9, qiskit 2.5.2,
pytket 2.18.1, Windows 11 x64, compactq 0.1.0 + native kernels 0.1.0) and are
reproducible with the commands shown. `gates / 2q / depth`.

Synthetic suite (`python -m compactq.bench`; same input QASM, same basis, fixed seeds):

| circuit | raw | Compact | qiskit L3 | 2q gain |
|---|---|---|---|---|
| ghz-5 | 5 / 4 / 5 | 5 / 4 / 5 | 5 / 4 / 5 | -0% |
| qft-3 | 18 / 6 / 14 | **14 / 6 / 11** | 14 / 6 / 11 | -0% |
| qft-4 | 34 / 12 / 22 | **25 / 12 / 17** | 25 / 12 / 17 | -0% |
| clifford-ladder-4 | 12 / 7 / 8 | 12 / 7 / 8 | 12 / 7 / 8 | -0% |
| clifford-ladder-5 | 15 / 8 / 13 | **13 / 8 / 11** | 13 / 8 / 11 | -0% |
| brickwork-4x4 | 44 / 12 / 16 | 40 / 12 / 17 | 40 / 12 / 16 | -0% |
| random-4q-40 | 40 / 18 / 24 | **37 / 18 / 21** | 40 / 22 / 27 | **-18%** |
| random-5q-60 | 60 / 22 / 37 | **49 / 22 / 30** | 62 / 28 / 37 | **-21%** |

Real circuits — QASMBench **small** suite (`python scripts/realbench.py`; unitary
cores only — circuits with classical control flow or mid-circuit measurement are out
of scope for a unitary optimizer and are skipped with a stated reason; 31 of 42 ran).
Every compactq row ≤ 8q is unitary-verified before being reported (0 inequivalent
outputs, refereed by `scripts/bench_json.py`); Qiskit L3 and pytket are unverified.
QASMBench assets are vendored in-tree (`third_party/QASMBench` — the small suite
plus the medium/large circuits tabulated below, pinned at pnnl/QASMBench `357b942`,
attribution in its LICENSE/NOTICE), so every number here is reproducible from a
fresh clone with no submodule step.

| circuit | Compact | qiskit L3 | pytket FullPeephole |
|---|---|---|---|
| adder_n10 | **110 / 57 / 95** | 137 / 65 / 99 | 165 / 61† / 118 |
| adder_n4 | **16 / 7 / 9** | 23 / 10 / 11 | 25 / 10 / 14 |
| basis_change_n3 | **34 / 10 / 22** | 49 / 10 / 28 | 79 / 10 / 50 |
| basis_test_n4 | 34 / 12 / 15 | 34 / 6¹ / 12 | 42 / 5† / 16 |
| basis_trotter_n4 | 773 / 240 / 352 | 794 / 179 / 361 | 969 / 159 / 419 |
| bell_n4 | **18 / 5 / 7** | 27 / 5 / 11 | 29 / 5 / 12 |
| cat_state_n4 | 4 / 3 / 4 | 4 / 3 / 4 | 6 / 3 / 6 |
| deutsch_n2 | 4 / 1 / 3 | 4 / 1 / 3 | 9 / 1 / 6 |
| dnn_n2 | **12 / 3 / 8** | 20 / 3 / 13 | 29 / 3 / 17 |
| dnn_n8 | **216 / 64 / 37** | 345 / 64 / 60 | 464 / 64 / 75 |
| error_correctiond3_n5 | **23 / 8 / 13** | 91 / 35 / 65 | 40 / 9 / 20 |
| fredkin_n3 | 19 / 8 / 11 | 19 / 8 / 11 | 22 / 8 / 13 |
| grover_n2 | 7 / 2 / 5 | 7 / 2 / 5 | 9 / 1† / 7 |
| hhl_n7 | **191 / 72 / 128** | 254 / 92 / 168 | 421 / 92 / 310 |
| hs4_n4 | 12 / 4 / 5 | 12 / 4 / 5 | 14 / 2† / 7 |
| ising_n10 | **166 / 49 / 29** | 260 / 90 / 46 | 370 / 90 / 58 |
| iswap_n2 | **7 / 2 / 5** | 8 / 2 / 6 | 9 / 1† / 7 |
| linearsolver_n3 | **11 / 4 / 9** | 17 / 4 / 12 | 25 / 4 / 19 |
| lpn_n5 | 7 / 2 / 4 | 7 / 2 / 4 | 15 / 2 / 8 |
| pea_n5 | **34 / 10 / 21** | 51 / 17 / 32 | 62 / 17 / 41 |
| qaoa_n6 | **114 / 36 / 47** | 166 / 36 / 63 | 197 / 36 / 83 |
| qec_en_n5 | 23 / 10 / 15 | 23 / 10 / 15 | 24 / 8† / 14 |
| qft_n4 | **20 / 6 / 10** | 34 / 12 / 20 | 39 / 12 / 23 |
| qrng_n4 | 4 / 0 / 1 | 4 / 0 / 1 | 12 / 0 / 3 |
| quantumwalks_n2 | **8 / 2 / 5** | 20 / 3 / 13 | 36 / 3 / 22 |
| sat_n7 | **125 / 52 / 71** | 158 / 60 / 86 | 193 / 60 / 107 |
| simon_n6 | **30 / 12 / 22** | 43 / 14 / 27 | 58 / 14 / 37 |
| teleportation_n3 | **5 / 2 / 4** | 6 / 2 / 4 | 12 / 2 / 8 |
| toffoli_n3 | **14 / 5 / 10** | 18 / 6 / 12 | 21 / 6 / 14 |
| variational_n4 | **29 / 8 / 14** | 44 / 8 / 18 | 51 / 8 / 24 |
| vqe_n4 | **25 / 9 / 11** | 46 / 9 / 18 | 75 / 9 / 23 |
| wstate_n3 | **18 / 6 / 12** | 22 / 6 / 15 | 30 / 6 / 21 |

† pytket 2.18.1's default-mode output fails unitary verification on this circuit
(fresh referee, `scripts/referee_pytket.py`: 6 of 30 ≤8q outputs inequivalent —
basis_test_n4 fid 0.50, grover_n2 0.50, hs4_n4 0.25, iswap_n2 0.50, qec_en_n5 0.25,
sat_n7 0.016; adder_n10 exceeds the 8q referee ceiling) — those counts are not
comparable. ¹ Qiskit L3 elides 2 SWAPs into `qc.layout` metadata here; reified they
cost 6 CX-equivalents — parity with compactq's 12.

Tallies vs Qiskit L3 (32 comparable rows from `python scripts/bench_json.py`, same
run): **2q win-or-tie 30/32** (11 wins / 19 ties / 2 losses), total-gate win **23 / tie
9 / lose 0**, depth win **22 / tie 9 / lose 1**. Zero compactq outputs flagged
INEQUIVALENT by the referee; Qiskit L3 is referee-flagged on basis_test_n4 (the
SWAP-elision above). The two 2q accounts not won: basis_test_n4 (12 vs 6 raw —
accounting parity once Qiskit's 2 elided SWAPs are reified) and basis_trotter_n4
(240 vs 179 — a pi/4-quantized Clifford+T ring whose optimization needs a
CliffordSimp-class fragment-resynthesis pass; on the roadmap).

Vs **pytket 2.18.1** (fresh referee, same protocol, `scripts/referee_pytket.py`):
default mode is inequivalent on **6 of 30** refereed circuits; of its 24 *valid*
outputs compactq wins the 2q count on 8 / ties 15 / loses 1 (basis_trotter_n4 —
the roadmap item above). In pytket's **safe** (`allow_swaps=False`) mode — now
itself inequivalent on 2 of 30 (basis_trotter_n4 fid 0.5, sat_n7 fid 0.016) — the
valid-output tally is **compactq 9 / tie 19 / lose 0**: zero valid 2q losses.

QASMBench **medium** (≤ 12 qubits; `python scripts/realbench.py --size medium
--max-qubits 12` — most medium files are ≥ 14q or carry classical control flow):

| circuit | Compact | qiskit L3 | pytket |
|---|---|---|---|
| sat_n11 | **486 / 212 / 355** | 599 / 252 / 403 | 713 / 250 / 507 |

(sat_n11 improved again this run: 212 2q vs the 252 of the 2026-09-08 run and a
252-vs-252 tie before that — the cross-pair/KAK pipeline keeps finding more.)

QASMBench **large** (≤ 32 qubits; `python scripts/realbench.py --size large
--max-qubits 32`; > 8q runs the no-verify path, soundness of which is fuzzed
against a Qiskit referee up to 10 qubits in the gauntlet):

| circuit | Compact | qiskit L3 | pytket |
|---|---|---|---|
| adder_n28 | **328 / 171 / 185** | 412 / 195 / 189 | 526 / 183 / 231 |
| bv_n30 | **55 / 18 / 20** | 79 / 18 / 20 | 197 / 18 / 24 |
| knn_n31 | **226 / 90 / 96** | 290 / 105 / 125 | 365 / 105 / 129 |
| qft_n29 | **805 / 370 / 85** | 1261 / 602 / 194 | 1699 / 806 / 197 |

(cc_n32 and vqe_uccsd_n28 skipped: classical control flow / non-standard QASM.
adder_n28 flipped from a 2026-09-08 loss (422/195 vs 412/195) to a clear win
(328/171 vs 412/195) with the current pass pipeline.)

MQT Bench — freshly generated algorithm-level circuits
(`python scripts/mqtbench_run.py`; mqt-bench, 44 circuits run, every tool's output
refereed at ≤8q; all four tools exact unless marked ‡). gates / 2q / depth:

| circuit | Compact | qiskit L3 | pytket FullPeephole | cirq CZTargetGateset |
|---|---|---|---|---|
| ghz_n4 | 4 / 3 / 4 | 4 / 3 / 4 | 6 / 3 / 6 | 11 / 3 / 7 |
| ghz_n6 | 6 / 5 / 6 | 6 / 5 / 6 | 8 / 5 / 8 | 17 / 5 / 11 |
| ghz_n8 | 8 / 7 / 8 | 8 / 7 / 8 | 10 / 7 / 10 | 23 / 7 / 15 |
| wstate_n4 | 13 / 6 / 8 | 13 / 6 / 8 | 25 / 6 / 14 | 21 / 6 / 11 |
| wstate_n6 | 21 / 10 / 12 | 21 / 10 / 12 | 41 / 10 / 20 | 35 / 10 / 17 |
| graphstate_n4 | 8 / 4 / 4 | 8 / 4 / 4 | 16 / 4 / 9 | 8 / 4 / 4 |
| graphstate_n6 | 12 / 6 / 4 | 12 / 6 / 4 | 24 / 6 / 12 | 12 / 6 / 4 |
| graphstate_n8 | 16 / 8 / 6 | 16 / 8 / 6 | 32 / 8 / 17 | 16 / 8 / 6 |
| qaoa_n4 | **20 / 4 / 10** | 23 / 8 / 15 | 50 / 8 / 33 | 25 / 8 / 17 |
| qaoa_n6 | **54 / 14 / 16** | 60 / 28 / 24 | 114 / 28 / 41 | 79 / 28 / 30 |
| qaoa_n8 | **100 / 28 / 36** | 108 / 56 / 54 | 182 / 56 / 66 | 163 / 56 / 70 |
| qft_n4 | **21 / 8 / 11** | 31 / 12 / 21 | 33 / 12 / 23 ‡ | 53 / 18 / 28 |
| qft_n6 | **44 / 18 / 17** | 71 / 30 / 35 | 73 / 30 / 37 ‡ | 110 / 39 / 45 |
| qft_n8 | **75 / 32 / 23** | 127 / 56 / 49 | 129 / 56 / 51 ‡ | 187 / 68 / 62 |
| qftentangled_n4 | **25 / 11 / 15** | 35 / 15 / 23 | 38 / 15 / 26 ‡ | 63 / 21 / 32 |
| qftentangled_n6 | **50 / 23 / 21** | 77 / 35 / 37 | 80 / 35 / 40 ‡ | 126 / 44 / 49 |
| vqe_real_amp_n4 | 25 / 9 / 11 | 25 / 9 / 11 | 73 / 9 / 26 | 33 / 9 / 15 |
| vqe_real_amp_n6 | 39 / 15 / 13 | 39 / 15 / 13 | 111 / 15 / 30 | 53 / 15 / 19 |
| vqe_real_amp_n8 | 53 / 21 / 15 | 53 / 21 / 15 | 149 / 21 / 34 | 73 / 21 / 23 |
| vqe_two_local_n4 | 34 / 18 / 17 | 34 / 18 / 17 | 82 / 18 / 29 | 59 / 18 / 28 |
| vqe_two_local_n6 | 69 / 45 / 25 | 69 / 45 / 25 | 141 / 45 / 37 | 130 / 45 / 44 |
| qnn_n4 | **11 / 3 / 5** | 20 / 3 / 8 | 39 / 3 / 15 | 13 / 3 / 7 |
| qnn_n6 | **17 / 5 / 7** | 29 / 5 / 9 | 59 / 5 / 19 | 21 / 5 / 11 |
| qpeexact_n4 | **18 / 5 / 10** | 20 / 7 / 13 | 25 / 7 / 15 ‡ | 30 / 10 / 21 |
| qpeexact_n6 | **51 / 16 / 20** | 64 / 27 / 37 | 76 / 27 / 42 ‡ | 96 / 33 / 49 |
| qpeinexact_n4 | **25 / 7 / 14** | 32 / 12 / 24 | 41 / 12 / 27 ‡ | 43 / 15 / 31 |
| qpeinexact_n6 | **52 / 17 / 22** | 72 / 30 / 44 | 87 / 30 / 47 ‡ | 103 / 36 / 56 |
| qwalk_n4 | **250 / 108 / 194** | 262 / 114 / 200 | 288 / 114 / 224 | 316 / 114 / 212 |
| qwalk_n6 | **1529 / 715 / 1261** | 1838 / 798 / 1382 | 2131 / 798 / 1653 | 2246 / 798 / 1416 |
| randomcircuit_n4 | **86 / 42 / 66** | 125 / 53 / 96 | 109 / 41 / 82 | 157 / 53 / 102 |
| randomcircuit_n6 | **180 / 78 / 104** | 229 / 95 / 134 | 272 / 94 / 151 | 269 / 95 / 146 |
| randomcircuit_n8 | **389 / 170 / 170** | 522 / 208 / 247 | 582 / 207 / 272 ‡ | 657 / 208 / 289 |
| cdkm_ripple_carry_adder_n4 | **27 / 15 / 25** | 34 / 17 / 26 | 43 / 16 / 30 | 48 / 17 / 32 |
| cdkm_ripple_carry_adder_n6 | **53 / 29 / 48** | 67 / 33 / 50 | 86 / 31 / 59 | 92 / 33 / 59 |
| draper_qft_adder_n4 | **15 / 5 / 11** | 20 / 8 / 17 | 25 / 8 / 23 | 23 / 8 / 17 |
| draper_qft_adder_n6 | **37 / 12 / 19** | 48 / 21 / 35 | 60 / 21 / 40 | 59 / 21 / 38 |
| grover_n4 | **97 / 44 / 77** | 135 / 52 / 95 | 158 / 52 / 109 | 143 / 52 / 95 |
| grover_n6 | **897 / 412 / 713** | 1098 / 456 / 818 | 1317 / 456 / 998 | 1287 / 456 / 810 |
| bv_n5 | 7 / 2 / 4 | 7 / 2 / 4 | 18 / 2 / 8 | 7 / 2 / 4 |
| bv_n7 | **4 / 3 / 4** | 10 / 3 / 5 | 24 / 3 / 9 | 10 / 3 / 5 |
| dj_n5 | 13 / 4 / 6 | 13 / 4 / 6 | 27 / 4 / 9 | 13 / 4 / 6 |
| dj_n7 | 19 / 6 / 8 | 19 / 6 / 8 | 39 / 6 / 11 | 19 / 6 / 8 |
| ae_n4 | **35 / 11 / 22** | 41 / 12 / 28 | 78 / 12 / 51 | 55 / 18 / 37 |
| ae_n6 | **66 / 24 / 36** | 92 / 30 / 55 | 145 / 30 / 91 | 125 / 42 / 62 |

Tallies (2-qubit count, invalid competitor outputs excluded): compactq wins
**25 / ties 19 / loses 0** vs Qiskit L3, **25 / 19 / 0** vs Cirq, and **14 / 19 / 1**
vs pytket — whose default-pipeline output is inequivalent on **10 of 44** circuits
(‡; the same SWAP-elision pitfall the QASMBench referee shows on pytket 2.18.1;
Cirq's optimizer was exact on all 44). Reproducibility note: cirq's QASM importer
and pytket's converter order qubits per register rather than by global wire index,
so the harness flattens every input onto a single q register first — without that,
both tools silently optimize the wrong unitary on multi-register circuits.

Honest reading: on the unitary QASMBench suite compactq wins or ties Qiskit L3 on
2-qubit count on every circuit except basis_trotter_n4 (240 vs 179 — its
pi/4-quantized iSWAP-ring structure needs deep template analysis; cross-pair
merging took compactq from 16 to 12 and it is a named roadmap item). At larger
scales the gap widens: on qft_n29 compactq cuts the 2q count a further 38% below
Qiskit L3 (370 vs 602) at less than half the depth, on ising_n10 it is the only
tool that finds non-trivial 2q reductions (49 vs 90), and **adder_n28 flipped from
the 2026-09-08 loss (422/195 vs 412/195) to a clear win (328/171 vs 412/195)**.
On latency, Qiskit's Rust core is still the fastest raw optimizer on mid-size
circuits, but for the same job — an optimized circuit you can actually trust —
compactq's proof-included time wins the majority of head-to-head circuits, and
pytket needs external verification on every output (its default mode is
inequivalent on 10 of 44 MQT and 6 of 30 refereed QASMBench circuits this run).

## CLI

Input policy: trailing measurements are dropped (the unitary core is optimized); mid-circuit measurement and reset are rejected with `UnsupportedCircuitError`; the Qiskit bridge rejects them outright (strip them first).  The CLI never crashes on width: circuits above the dense proof limit are routed to randomized verification, and every run reports its proof status

The console command is `compactq` on every platform — deliberately *not*
`compact`, which collides with Windows' built-in `compact.exe` file-compression
tool. `python -m compactq` always works too.

```bash
compactq in.qasm -o out.qasm        # optimize an OpenQASM 2.0 file
compactq in.qasm --stats            # print gate-count/depth deltas
compactq in.qasm --approx 0.99      # bounded-fidelity approximate mode
compactq in.qasm --objective depth  # depth-first: depth may never grow
python -m compactq.bench            # full benchmark vs Qiskit (if installed)
python -m compactq.bench --quick    # 3-circuit CI guard
compact-bench --out results.md     # installed console script
```

## API

```python
import compactq
compactq.optimize(circuit, verify=True)        # Circuit -> Circuit (smaller, verified)
compactq.optimize_deep(circuit)                # KAK cascade, best for dense circuits
compactq.optimize_search(circuit)              # multi-pipeline search, keeps the best
# objectives — which metric may never grow (all exact; lexicographic orders):
compactq.optimize(c, objective="2q")           #   default: (2q, gates, depth)
compactq.optimize(c, objective="depth")        #   depth-first acceptance
compactq.optimize(c, objective="gate_count")   #   total-gates-first
compactq.optimize(c, objective="latency")      #   depth alias
compactq.optimize_search(c, objective="weighted")  # min 1.0*2q + 0.1*depth + 0.02*gates
compactq.optimize_for(c, target)               # hardware-error weighting (compactq.target)
compactq.verify(original, optimized)           # independent evidence: {'equivalent', 'tier', 'method', ...}
compactq.approximate(circuit, min_fidelity=0.99)  # trade bounded fidelity for fewer 2q gates
from compactq.target import Target, optimize_for
t = Target(cx_fidelity={(0, 1): 0.999, (1, 0): 0.98})
best, est_fid = optimize_for(circuit, t)      # least-noisy circuit for YOUR machine
from compactq.target import approximate_for_target
best, est_fid = approximate_for_target(circuit, t)  # greedy fidelity-budget allocation
compactq.to_qasm(circuit) / compactq.from_qasm(t)  # OpenQASM 2.0 round-trip
compactq.from_qasm3(t)                         # OpenQASM 3 import (common subset)
compactq.to_qasm3(circuit)                     # OpenQASM 3.0 export
compactq.benchmarks.qft / ghz / brickwork / random_circuit / clifford_ladder
compactq.param / compactq.structure_optimize / compactq.bind   # symbolic angles: optimize the structure once, bind later
compactq.optimize_large(circ)   # verified optimization to ~30q (randomized K-state proof; numpy)
```

`from_qasm` understands the full extended-qelib1 set (u1/u2/u3/u, sx, sxdg, cy, cz,
swap, cswap, ccx, crz, cu1/cp, rzz, rxx, id), multiple registers with global wire
numbering, register-wide operands (`h q;`) and user `gate` definitions.
`from_qiskit` accepts any `QuantumCircuit`: gates inside compactq's IR pass through
untouched, everything else (`mcphase`, `ccx`, `ecr`, `iswap`, `cu`, `rxx`, ...)
is decomposed at the boundary into the supported basis via Qiskit's own
equivalence library — loss-free, and the result gets optimized instead of
crashing.

Hardware-aware passes (`compactq.hardware`):

```python
from compactq.hardware import route, route_aware, flip_cx, translate_1q_to_rz_sx_x
routed, final_map = route(circuit, coupling=[{0,1},{1,2}])
routed, final_map = route_aware(circuit, coupling, target=t, restore=True)
```

`route` inserts SWAPs along shortest paths and returns the final
logical-on-physical mapping (apply it to your measurements); `route_aware`
adds SABRE-style re-routing with error-weighted look-ahead and optional
permutation restore. Optional bridges:
`compactq.qiskit_bridge` (`compactq_pass`, `to_qiskit`, `from_qiskit`) and
`compactq.cirq_bridge` (`to_cirq`).

Native hardware bases: `Target(native_2q="ecr" | "cz" | "iswap")` with `optimize_for`, or `compactq --native ecr` on the CLI (verified gate-count parity with Qiskit's own basis decomposer).

Supported gates: `h x y z s sdg t tdg rx ry rz p sx sxdg u/u3 cx cz swap cp`,
plus `mcx`/`mcp` (expanded on use). Multi-controlled and exotic gates arriving
via `from_qiskit` are boundary-decomposed automatically.

## Project history

- **v0.1.0 — initial public release (2026-09-17).** The verified
  optimizer (peephole / commutation / CP / pure-Python KAK /
  Clifford-tableau / phase-polynomial passes, multi-pipeline search,
  exact + approximate modes, hardware targets, routing, QASM2/3 IO,
  Qiskit bridge + plugin, optional Rust kernels) plus the open
  error-suppression stack (twirling, benefit-gated DD, MLE mitigation,
  ZNE, CDR, simulators, device metrics, execution adapters). See
  [CHANGELOG.md](CHANGELOG.md) for the full inventory.

  The public version series starts at 0.1.0. An earlier private
  prototyping sprint (local iterations, 2026-09-10 → 2026-09-16)
  produced the engine; those internal numbers are retired and the git
  history keeps the full trail.

## Roadmap

**Research order** (post-release review, 2026-09-17): 1. Trotter/Hamiltonian
optimization → 2. scalable verification → 3. BQSKit benchmark
(`scripts/bqskit_bench.py` ships; runs once BQSKit is installed) →
4. hardware-mapped benchmark → 5. random-SU(4) synthesis study (10,000
blocks, KAK vs Qiskit/TKET/BQSKit) → 6. large randomized circuits →
7. calibration-aware optimization.  Simulators and PEC come after — they
must not distract from the compilation thesis.

- **Primary research target — Trotter / Hamiltonian simulation.**
  `basis_trotter_n4` (240 vs Qiskit's 179 2q) is the flagship measured
  loss, and large structured phase/Trotter rings are exactly where
  external benchmarking (Benchpress, *Nature Computational Science*
  2025) reports the biggest competitor gains.  Matching them requires
  layer-boundary content-permutation search with swap-network payoff at
  the output (the Clifford+T normal-form pass in
  `compactq/cliffordt.py` covers only the pi/4-quantized ring; the real
  trotter file contains arbitrary-angle PhasedISWAPs and needs a
  different attack).  qiskit's 179 and pytket's 159 are both
  *permutation-elided* (unverified outputs; their refereed fidelity on
  this circuit was 0.25/0.5) — our 240 is exact.
- **Scalable verification.** Evolve `optimization → full-unitary proof →
  return` into local proof certificates + compositional equivalence, so
  the verified regime extends past 8 qubits structurally (randomized
  K-state verification covers ~30 qubits today).
- **Benchmark suite v1.** Add a measured **BQSKit** column and the four
  MQT Bench abstraction levels (algorithmic / target-independent /
  native-gate / hardware-mapped) to the harness; publish every raw
  artifact under `results/` with environment metadata.
- **Calibration-aware compilation** (generalizes the fixed `weighted`
  weights and the depth-as-latency proxy): cost = α·2q + β·depth + γ·1q +
  δ·hardware_error + ε·duration, with the coefficients coming from the
  `Target` calibration rather than being fixed constants.  Note that
  *compilation* latency and *execution* latency are different objectives;
  today's `latency` objective is the execution-depth proxy.
- Rust parity-network BFS kernel: packed
  u64 wire-mask states; exactness fuzzed against the phase-polynomial
  reference (0 failures); parity_pass 302ms -> 25ms on parity-heavy
  shapes; finds windows the Python budget quirk misses
- port the pass-engine hot path (peephole/KAK driving loops) to Rust to
  close the remaining raw-latency gap vs Qiskit on large circuits
- MQT Bench in CI (currently a local script: `scripts/mqtbench_run.py`)

## Development

```bash
git clone https://github.com/Q-PROOF/Compact && cd Compact
python tests/run_tests.py     # zero-dependency test suite
python scripts/gauntlet.py    # 1,018-check end-to-end gauntlet (needs qiskit)
pip install -e .[bench]       # optional: qiskit for the comparison column
python scripts/realbench.py --size small            # QASMBench vs Qiskit/pytket
python scripts/realbench.py --size large --max-qubits 32
pip install mqt-bench cirq ply   # extras for the fourth suite, then:
python scripts/mqtbench_run.py   # MQT Bench: compactq vs Qiskit/pytket/Cirq
```

Building the optional Rust wheel:

```bash
cd native && pip install maturin && maturin build --release -o dist
pip install dist/compactq_native-*.whl
```

MIT licensed. Contributions welcome — every PR must keep the property tests green.
