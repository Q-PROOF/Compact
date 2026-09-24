# CORRECTNESS.md — break our prover, get credited

Compact's core claim is: **every returned circuit is unitarily
equivalent to the input (up to global phase), and every equivalence
verdict is backed by the proof tiers defined in
[VERIFICATION.md](VERIFICATION.md).**  This document asks you to break
that claim — and commits to making the process fast and public.

## The bounty invitation

If you can construct **any** circuit where

1. `optimize` / `optimize_search` / `optimize_large` / `optimize_for`
   returns an output that is NOT unitarily equivalent to the input
   (up to global phase), or
2. `compactq.verify(a, b)` answers `equivalent: true` for circuits
   that differ, or answers `equivalent: false` for circuits that are
   equal (false inequivalence is also a bug — it must decline instead,
   see VERIFICATION.md), or
3. a certificate passes `compactq-check` for circuits that differ, or
4. the documented tier/boundary in VERIFICATION.md does not match what
   the code actually does,

then file an issue with the **smallest QASM reproducing it**:

**https://github.com/Q-PROOF/Compact/issues** — label `correctness`.

Confirmed correctness bugs are treated as the highest priority in the
repository, fixed before any feature work, credited in the CHANGELOG
("found via adversarial testing, by @you"), and, if you opt in, added
to the permanent regression suites so your case can never regress.

## Triage process

- **Acknowledge:** within a few days of filing.
- **Reproduce or reject:** we run your QASM through the current release
  and reply with either a confirmed reproducer or the exact
  transcript showing the claimed behavior does not reproduce.
- **Fix + regression test:** every accepted report lands with a named
  test in `tests/` (your choice of name credit) before the fix release.
- **Wrong reports are fine.**  A report that doesn't reproduce costs
  you nothing; sloppy handling of a *correct* report costs us
  everything.

## What we red-team ourselves (and what it found)

Before each release we attack the provers with adversarial suites:
`tests/test_adversarial_v024.py` (and the historical cases in
`tests/run_tests.py::test_adversarial`) throw empty circuits, identity
chains, global-phase traps, degenerate angles (0, 2π, 1e-18, 100π),
tier-straddling Clifford+T, node-budget boundaries, cross-prover
verdict-agreement sweeps, and malformed inputs at the stack.

Found by adversarial testing so far:

- **v0.2.4** — `cx` on a repeated wire (`cx q[0], q[0]`) was silently
  accepted by the IR; a gate with no well-defined unitary must never
  enter a proven circuit.  Now rejected at construction.
  *(self-reported, fixed)*
- **v0.2.4** — qelib1's `ch` (controlled-H) failed QASM import with
  `unsupported 2q gate`, against the extended-import contract.  Now
  imported via its exact 8-gate decomposition (verified against the
  dense matrix, both wire orders).  *(self-reported, fixed)*
- **v0.2.4** — the decision-diagram prover returned a false
  INEQUIVALENT verdict (tier 2!) comparing a ≥7-control MCX circuit
  against its own copy: near-cancelling blocks amplified under
  first-nonzero normalization and corrupted the diagram.  Fixed by
  max-weight (QMDD-standard) normalization plus a self-consistency
  norm guard that converts any residual corruption into a loud
  decline.  *(self-reported, fixed — see `tests/test_adversarial_
  v024.py::test_dd_self_comparison_mcx`)*
- **v0.2.4** — the independent checker's dense-unitary kernel rebuilt
  the source row from the COLUMN's other bits and skipped recomputing
  cells the gate did not mix, so unitary-witness verdicts on any
  circuit with off-diagonal structure were wrong (fidelity 0.68
  reported for a valid QFT-4 certificate); its Clifford canonical form
  also crashed on tuple generators.  Both were pre-existing (present in
  the pristine v0.2.3 checkout — the suite was not CI-registered).
  Fixed: correct row-based kernel, full-cell recomputation, list
  generators; `tests/test_checker_independence.py` is now in CI.
  *(self-reported, fixed)*
- **v0.2.4** — `verify()` dispatched to randomized statevector
  verification with no width ceiling; a 33q circuit would attempt a
  2³³-state allocation and hang instead of declining.  Now capped at
  the documented 30q memory ceiling (`verify_large.RAND_MAX_QUBITS`).
  *(self-reported, fixed)*

## Scope notes (so we agree on the rules)

- Approximate mode (`compactq.approximate`, `optimize_for` with a
  fidelity target) trades fidelity **by explicit, measured budget** —
  reporting that a result is approximate is not a bug; a result that
  misses its *measured bound* is.
- Equivalence is up to global phase (the physics standard); reporting
  `optimize` "changed" a circuit by a global phase alone is not a bug.
- Suppression/DD-insertion and mitigation layers are
  simulator-evidence-based, not exactness-proven — their claims live in
  `docs/suppression.md` and the `results/` artifacts, not under this
  bounty's exactness contract.
