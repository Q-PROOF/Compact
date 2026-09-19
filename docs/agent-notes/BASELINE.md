# BASELINE — recorded before agent implementation work

- commit: 8934bfe → 371ce35 line (v0.2.0 release), pre-work HEAD to be stamped
  by the first agent commit on top
- python: 3.11.9 (Windows 11 x64, pyenv-win)
- compactq: 0.2.0 (PyPI-published; working tree at v0.2.0 + provenance/CI commits)

## Gates (all green at baseline)

- `python tests/run_tests.py` — PASS (76 functions)
- `python tests/test_objectives.py` — PASS (4)
- `python tests/test_verify_api.py` — PASS (5)
- `python tests/test_trust_layer.py` — PASS (5+2)
- `python scripts/gauntlet.py` — 1018/1018, 0 inequivalent
- `python scripts/version_lint.py` — OK
- `python scripts/check_provenance.py` — OK (all 8 artifacts carry provenance)
- `python -m compactq.bench --quick` — zero regressions
- `python scripts/bench_gate.py` — PASS (32 circuits; depth allowance now
  cross-platform scaled: max(2, 5% of baseline depth))

## Known reference numbers (from committed artifacts)

- QASMBench small (43 rows, referee: Qiskit Operator): compactq outputs
  exact-unitary on all ≤8q rows; 0 inequivalent anywhere
- independent referees: PyZX 27/27 verified (2 parser-quirk errors
  recorded); MQT QCEC 11 OK + 2 tolerance-disagreements (Operator
  fidelity 1.0 / 0.999999999972) + 11 non-unitary skips
- BQSKit head-to-head (32 circuits): compact 2q win 10 / tie 19 /
  lose 3; geomean bqskit/compact 1.16; 0 invalid bqskit outputs
- scalability: 112 records, 4→256Q; ceilings (strict practical class):
  correctness/stability 256Q, runtime 256Q, memory 256Q
- repro harness (17 independent circuits, u3+cx lowered): compact mean
  2q cut 10.2% (median 21 ms), qiskit L3 14.6% (4 ms), pytket 25.0%,
  cirq 2.8%; compact vs qiskit 2q W/T/L 2/12/1
- known measured loss: grover_5 386 vs 144 2q (wide-MCX bridge
  expansion cost — roadmap item)
- external audit (independent agent, 2026-09-18): 0 wrong compactq
  outputs on 360 verify=False runs + refereed suite; CX-lowered
  head-to-head 0/34/3 (ties dominate); suppression gains
  noise-model-dependent (×1.04–1.09 in external decoherence models)

Every later phase is measured against these numbers. The 2q count is the
hard invariant in `scripts/bench_gate.py`; gates/depth carry documented
cross-platform drift allowances.
