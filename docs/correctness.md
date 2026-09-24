# Correctness model

> **This document moved.**  The authoritative, tier-by-tier
> methodology is **[VERIFICATION.md](../VERIFICATION.md)** (root of the
> repository, linked from the README) — it defines every proof tier,
> the exact boundary of each prover's coverage (each backed by a named
> test), and the loud-decline contract.  The adversarial-testing
> invitation lives in **[CORRECTNESS.md](../CORRECTNESS.md)**.

Short form:

- Every rewrite is **exact**: the unitary is preserved up to global phase.
- Proofs are attached, not assumed: `compactq.verify(a, b)` returns the
  verdict plus the evidence tier (0–4), and `optimize_search()` /
  `optimize_large()` return your input unchanged whenever the proof
  layer declines.
- Approximate mode (`compactq.approximate`) is the one deliberate
  exception: each replaced block carries a **measured** fidelity
  guarantee, and total infidelity is bounded by (#approximated blocks)
  x (1 - min_fidelity).

## Test suite

```bash
python tests/run_tests.py                # zero dependencies required
python tests/test_adversarial_v024.py    # adversarial red-team probes
python tests/test_tier_boundaries.py     # VERIFICATION.md boundary pins
python tests/test_transpiler_passes.py   # qiskit/pytket integration (extras)
```

The suite: unit tests, regression tests for classic optimizer bugs,
oracle comparisons (qiskit, pytket), randomized property tests, the
adversarial red-team probes, and the tier-boundary pins that keep
VERIFICATION.md and the code from drifting apart.
