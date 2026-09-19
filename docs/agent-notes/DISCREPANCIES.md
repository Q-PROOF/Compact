# DISCREPANCIES (repo vs ForSure-Implementation.md)

- §Baseline says "v0.1.0 released 2026-09-17, 0 stars": the live repo is at
  v0.2.0 (0.1.0→0.1.8 shipped 2026-09-17..19) with artifacts on PyPI and
  GitHub Releases. The plan's version ladder is shifted accordingly
  (certificates land in 0.2.1 rather than being "0.2.0").
- §T0.3 basis sets: implemented with `{cz, rz, sx, x}` and the ECR variant
  via `compactq.target` native_2q ("cz"/"ecr") — matches the intent.
- §6.1 certificate example field `"composition": "tensor_local_v1"`: v1
  whole-circuit certificates use `"composition": "whole_circuit_v1"`
  (regions arrive with T1.2). Spec updated to match the implementation.
- §12.4 publish order (`compactq-check` first): adopted.
