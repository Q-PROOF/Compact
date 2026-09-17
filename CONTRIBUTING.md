# Contributing to Q-PROOF Compact

Thanks for your interest! The project's core promise is simple: **every
optimized circuit ships with a correctness argument**. Any contribution
that weakens that promise will be declined, no matter how fast or how
much it reduces gates.

## Getting started

```bash
git clone https://github.com/Q-PROOF/Compact && cd Compact
python tests/run_tests.py        # zero-dependency test suite
pip install -e .[bench]          # optional: qiskit for oracle tests
python scripts/version_lint.py   # docs/metadata consistency
python scripts/gauntlet.py       # end-to-end gauntlet (needs qiskit)
```

## The rules any PR must keep

1. **Never ship unproven.** Every rewrite must be exact up to global
   phase and covered by a proof path (dense unitary, tableau, per-block
   measured fidelity, or randomized K-state verification). If your pass
   can't prove it, don't accept it — return the input unchanged.
2. **Run the gates.** `tests/run_tests.py` green, `scripts/gauntlet.py`
   green, and benchmarks re-run if you touched a pass that can change
   gate counts. Benchmark claims must come from measured runs, never
   hand-edited numbers.
3. **Zero-dependency core.** `compactq/` must import without numpy,
   qiskit, or any third-party package. Optional accelerators
   (`compactq_native`, numpy paths) must degrade gracefully and are
   validated separately.
4. **Scope honesty.** The optimizer transforms unitary circuits.
   Non-unitary operations (mid-circuit measurement, reset, classical
   control flow) are rejected with `UnsupportedCircuitError` — never
   silently dropped. Don't add silent coercion anywhere.

## What's worth contributing

- Algorithmic passes (see README roadmap: layer-boundary permutation
  search, T-aware CX hoisting, Rust pass-engine loops)
- Benchmark/QASMBench/MQT reproducibility improvements
- Documentation, examples, and external-oracle test coverage
- Native-kernel coverage for more pass hot paths

## Reporting bugs

Open an issue with: the exact circuit (QASM preferred), the entry point
used, the version (`compactq.__version__`), and — critically — **the
proof status reported**. A circuit that comes out inequivalent is a
release blocker; report it even if you already found the workaround.

## Security

See `SECURITY.md` for reporting policy.

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
By participating, you agree to uphold its standards.  Report violations
to **qproof.ai@gmail.com**.
