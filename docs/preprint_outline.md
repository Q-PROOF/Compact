# Preprint outline — "Compact: a zero-dependency, proof-carrying circuit
# optimizer with tiered verification"

Working notes toward an arXiv preprint (quant-ph).  Every claim below
maps to a committed artifact in `results/` — the paper must not contain
a number that cannot be regenerated from a producer script.

## Title

Compact: verified quantum circuit optimization with tiered,
machine-checkable evidence

## Abstract (skeleton)

- Problem: production optimizers silently transform circuits with no
  correctness evidence; incorrect optimized circuits have been observed
  in widely used default pipelines.
- Contribution: a zero-dependency optimizer whose every output carries
  tiered, machine-checkable evidence of unitary equivalence up to global
  phase, plus an independent verification CLI and certificate format.
- Results summary (regenerate at submission time): correctness on
  refereed benchmark sets; independent cross-refereeing (Qiskit
  Operator, PyZX ZX-tensors, MQT QCEC decision diagrams); stability
  frontier to 256 qubits; head-to-head against Qiskit L3, pytket, Cirq
  and BQSKit.

## 1. Introduction

- Correctness gap in production quantum compilers; cite the observed
  inequivalent default-pipeline outputs and the SWAP-elision pitfall.
- Thesis: optimization and evidence are one artifact.

## 2. Design

- Tiered verification ladder T0-T5 (unverified, verified-randomized,
  exact-proven numerical, exact-proven algebraic, compositional, formal).
- Proof-carrying pipeline: verify-before-return; decline on doubt.
- Certificate format: hashes, metrics, evidence, reproducibility.

## 3. The optimizer

- Pass inventory: peephole, commutation, CP/phase-polynomial, KAK/Weyl,
  Clifford tableaux, search; multi-objective acceptance orders.
- Zero-dependency core; optional native kernels.

## 4. Methodology

- Same-input, same-gate-set harness; independent referees (Qiskit
  Operator, PyZX, MQT QCEC); provenance blocks; reproducibility gate.
- Honesty rules: losses reported; measurement conditions documented.

## 5. Results

- QASMBench small/medium/large; MQT Bench; independent-circuit
  reproduction; BQSKit head-to-head; scalability frontier to 256Q.
- Tables regenerate from `results/` artifacts.

## 6. Limitations

- Numerical certificates are not formal proofs; equivalence-proven
  ceiling is structure-dependent at large widths; runtime frontier on
  dense structured workloads; simulator-only suppression evidence.

## 7. Related work

- Qiskit transpiler, TKET, BQSKit, PyZX, MQT QCEC/Bench, VOQC/VyZX,
  Quasar (equality saturation), FT compilation (Litinski, gridsynth,
  PBC), FTCircuitBench.

## 8. Conclusion and outlook

- T5 formal extraction for the algebraic fragment; hardware-validated
  suppression; 3Q/4Q synthesis.
