# Correctness model

- Every rewrite is **exact**: the unitary is preserved up to global phase.
- For circuits within the dense proof limit (8 qubits with the optional native wheel, 6 without),
  `optimize_search()` **proves** equivalence via full-unitary comparison
  before returning; on any numerical doubt it returns your circuit unchanged.
- Clifford blocks are proven with stabilizer tableaux — an exact proof valid
  at **any** qubit count, independent of the unitary limit.
- The KAK/Weyl re-synthesizer is cross-validated against Qiskit's
  implementation (Weyl coordinates agree to ~1e-15 on randomized SU(4)s) and
  every re-synthesized block is re-verified against its own 4x4 unitary
  before it can replace anything.
- The CX update rules of the tableau engine are verified exhaustively against
  matrix conjugation over all phase states and both CX directions.
- Approximate mode (`compactq.approximate`) is the one deliberate exception:
  each replaced block carries a **measured** fidelity guarantee, and total
  infidelity is bounded by (#approximated blocks) x (1 - min_fidelity).

## Test suite

```bash
python tests/run_tests.py     # zero dependencies required
```

The suite (50+ checks and growing): unit tests, regression tests for classic optimizer bugs, oracle
comparisons (qiskit, pytket), and randomized property tests.
