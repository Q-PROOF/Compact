# DECISIONS

| id | decision | rationale | date |
|---|---|---|---|
| D1 | Agent commits directly to `main` (no PR branches) | the operator's established pipeline is commit → full gate chain → push to main; §2.3 branch/PR flow is adopted for external contributors instead | 2026-09-19 |
| D2 | Prover = Lean 4 preferred; **deferred** to a dedicated phase | no Lean toolchain in the working environment; §T1.3 fallback allows shipping the enumerated rule library + unitary-verified tier first; recorded as `formal/RULES.md` + CUTS | 2026-09-19 |
| D3 | `compactq-check` v1 witnesses are **whole-circuit** (unitary ≤8q, stabilizer any-width, phase-polynomial any-width); per-rewrite regions deferred to T1.2 | per-rewrite region instrumentation across every pass is a multi-week instrumentation effort; whole-circuit witnesses already deliver M1 (independent re-checkability) soundly; region decomposition follows in the next phase | 2026-09-19 |
| D4 | Certificates are opt-in via `compactq.cert.optimize_with_certificate()`; `optimize_search` signature unchanged | backward compatibility per §T1.1.2; the result-object kwarg form lands with region support | 2026-09-19 |
| D5 | CI matrix adds windows/macos for the zero-dep test job; qiskit regression gate stays ubuntu-only | cross-platform 1q tie-break drift is documented and allowance-covered (bench_gate); full matrix on the qiskit job would triple CI minutes for no new signal | 2026-09-19 |
| D6 | Quasar/Quarl dispositions recorded as NOT_RUNNABLE, never estimated | §T0.2 hard rule | 2026-09-19 |
