# compactq gym for Qiskit Benchpress

A [Benchpress](https://github.com/qiskit/benchpress) gym for compactq —
the upstreamable form of our Benchpress integration: one new package
(`compactq_gym/`) plus one `elif` branch in each of Benchpress's two
per-gym dispatch files.  These are exactly the files an upstream PR
would add/change.

## What it measures

The Benchpress **abstract-transpilation workout** on the QASMBench
suites at the all-to-all topology, through Benchpress's own pytest
harness and `pytest-benchmark` records — same circuits, same loader,
same record schema as the qiskit gym, so rows are directly comparable.
The timed region for compactq is exactly `optimize_search` — its
whole-circuit equivalence proof is included in the time, because
compactq never returns unproven output.

## Run it

```bash
pip install pytest pytest-benchmark
python scripts/benchpress_run.py --gym compactq --sizes small,medium
python scripts/benchpress_run.py --gym qiskit   --sizes small,medium   # paired baseline
```

Records land in `results/benchpress/` (native pytest-benchmark JSON,
committed unchanged).  `BENCHPRESS_HOME` can point at an existing
clone; otherwise a pinned github tarball is downloaded to
`third_party/benchpress` (gitignored) with host-allowlisted fetch and
traversal-checked extraction.

## Files and what an upstream PR would do with them

| file | upstream action |
|---|---|
| `compactq_gym/**` | add as `benchpress/benchpress/compactq_gym/` |
| `patches/apply_patches.py` | replaced by the two `elif` insertions it performs in `utilities/io/circuit_output.py` and `utilities/validation/validation.py` (applied idempotently here because we cannot patch upstream from a release) |

Scope note: the compactq gym claims the **all-to-all** abstract
topology only (logical optimization; routing is a separate compactq
pipeline).  Equivalence of every output is proven in-product before
compactq returns — the Benchpress validator checks the output
contract, not correctness.
