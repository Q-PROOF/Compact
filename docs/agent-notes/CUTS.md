# CUTS — implemented subsets and deferred scope

Per §0.2: infeasible tasks are implemented as the closest sound subset
and recorded here. Nothing on this list is silently dropped.

| task | cut / subset | why | follow-up |
|---|---|---|---|
| T1.2.2 per-rewrite region emission in every pass | v1 certificates certify the whole circuit (unitary / stabilizer / phase-polynomial witnesses); no `rewrite_trace` yet | instrumenting every pass with region hooks is a multi-week change across the entire pass inventory | T1.2 (next phase) |
| T1.3 Lean formalization of ≥20 rules | `formal/RULES.md` enumerates the 20-rule set with statements; Lean project deferred (D2) | no Lean toolchain in environment | dedicated phase; CI job `formal` when added |
| T1.1.2 `optimize_search(certificate=True)` kwarg | opt-in `compactq.cert.optimize_with_certificate(circ)` returns the result object; default `optimize_search` signature untouched | backward compat per §T1.1.2; kwarg form lands with region support (D4) | T1.2 |
| T0.2 Quarl comparator | NOT_RUNNABLE: no PyPI distribution (`quarl` not found) | — | revisit if released |
| T0.2 Quasar comparator | NOT_RUNNABLE: PLDI'26 research artifact (Zenodo), not a pip-installable tool in this environment | §T0.2 forbids estimating its numbers | revisit with artifact build |
| T0.5 Python 3.9 in CI matrix | matrix ships 3.10–3.13 × ubuntu/windows/macos; 3.9 added after a dedicated compat pass | `requires-python >=3.9` is declared but 3.9 was never exercised in this environment; adding an untested red matrix row would break the gate honestly but block all merges | compat pass, then extend matrix |
| Phase 2 (e-graph, rule generation, Rust engine, 3Q/4Q resynthesis) | roadmap ladder (0.3.0) | multi-week engineering campaigns; §0.2 subset rule | per §8 task groups |
| Phase 3 (T-count passes, PPM/Litinski, resource-estimator export) | roadmap ladder (0.4.0); `compactq.resources` T-count/T-depth + Clifford+T rebase already ship as the seed | §8 scope discipline: verified front-end first | per §8 task groups |
| Phase 4 (real-hardware suppression evidence) | harness design committed in README honesty notes; phase marked **BLOCKED: needs hardware credentials** per §9 hard rule | no IBM/Braket credentials in environment | the moment credentials exist |
| Phase 5 referee service/leaderboard | `compactq verify` CLI + certificate format ship now (the substrate); public leaderboard service deferred | service hosting out of scope for a code session | T5.2 |
| arXiv preprint | `docs/preprint_outline.md` committed (full skeleton + claims-to-evidence map) | writing/submission is the operator's call | T5.3 |
