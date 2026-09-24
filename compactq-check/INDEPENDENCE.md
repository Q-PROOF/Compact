# Implementation-independence audit — compactq-check vs compactq

**Audited:** 2026-09-23 (v0.2.4 release cycle) · **Result: INDEPENDENT**

## Why this matters

`compactq-check` exists so that a Compact certificate can be verified by
code that did not produce it.  If the checker shared a code path with
`compactq` — the same unitary engine, the same Clifford tableau update
rules, the same decision-diagram implementation — a single bug would
produce matching (wrong) answers on both sides and the check would be
worthless.  The audit question is: does the checker share ANY code path
with the `compactq` package?

## How independence is enforced

1. **No imports.**  `compactq-check/src/compactq_check/` imports only
   the Python standard library and its own modules (`checker.py`,
   `ops.py`, `qasm.py`, `cli.py`).  A grep for cross-package imports is
   part of the release checklist:

   ```
   grep -rn "import compactq\|from compactq" compactq-check/src/
   # -> no matches
   ```

2. **Runtime proof.**  `tests/test_checker_independence.py::test_
   checker_never_imports_compactq` purges every `compactq*` module from
   `sys.modules`, imports `compactq_check.checker`, and asserts nothing
   `compactq`-prefixed reappears.  Runs in CI on every push.

3. **Separate packaging.**  The checker ships as its own package
   (`compactq-check/pyproject.toml`), installable without `compactq`
   present.  Its QASM parser (`qasm.py`), unitary builder, Clifford
   canonicalization, phase-polynomial extraction, and decision-diagram
   overlap check are written fresh against the certificate
   specification (`docs/CERTIFICATE_SPEC.md`), not against compactq's
   source.

## What IS shared (deliberately)

- **The specification, not the code.**  Both sides implement
  certificate schema v1, the same witness kinds and limits
  (`unitary` ≤ 8q, `dd` ≤ 32q, `stabilizer`/`phase_polynomial` any
  width) and the same decision rule (Hilbert-Schmidt fidelity
  |Tr(A⁺B)|/d > 1 − 1e-7, equivalence up to global phase).  This
  agreement is the point: two implementations of one spec.
- **Tolerances.**  A shared tolerance is a spec-level choice; it makes
  verdicts comparable.  It is cross-checked from the outside by the
  referee layer (Qiskit `Operator`, MQT QCEC), which shares no code
  with either side.

## Residual risk, stated honestly

Shared-spec/different-implementation catches **implementation** bugs
(one side wrong, the other right → INVALID, investigation triggered).
It cannot catch **specification** bugs (both sides agree on a wrong
rule).  Spec-level claims are therefore additionally audited against
external referees — qiskit's Operator as the dense oracle and QCEC as
an independent decision-diagram-class prover — see `results/
qcec_referee.json` and `scripts/referee_qcec.py`.

## Change policy

Any change to `compactq`'s provers does **not** require a matching
checker change (and vice versa) unless `CERT_VERSION` is bumped — in
which case the checker must be updated in the same release, and the
release notes must say so.
