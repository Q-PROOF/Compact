# Pull Request

## What does this PR do?

<!-- Describe the change: new pass, bug fix, benchmark, docs, etc. -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New optimization pass
- [ ] New suppression/mitigation feature
- [ ] Benchmark / performance improvement
- [ ] Documentation update
- [ ] Test / CI improvement

## Checklist

- [ ] `python tests/run_tests.py` passes (zero dependencies)
- [ ] Every rewrite is exact up to global phase (proof-net verified)
- [ ] `compactq/` core imports without numpy/qiskit/third-party
- [ ] Benchmark claims come from measured runs, not hand-edited
- [ ] Version-consistency lint passes (`python scripts/version_lint.py`)
- [ ] If a new pass: property-tested (500+ random-circuit trials)

## Proof status

<!-- For optimization passes: what proof covers the output?
     exact-unitary / tableau / randomized-K / measured-fidelity -->

## Additional context

<!-- Any relevant benchmarks, circuit examples, or references -->
