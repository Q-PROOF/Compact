# Scalability gauntlet (generated - do not edit)

Generated 2026-09-18T16:30:27.137934+00:00 | commit `387f3576d010` | compactq 0.1.3 | referee limit 10q (independent Qiskit Operator); Clifford workloads are algebraically proven at any width.

| width | tests | crashes | timeouts | incorrect | policy viol | nondet | verified | median opt ms | peak MB | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 78 | 0.9 | PASS |
| 6Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 124 | 2.0 | PASS |
| 8Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 94 | 0.3 | PASS |
| 10Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 143 | 0.4 | PASS |
| 12Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 143 | 0.5 | PASS |
| 16Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 205 | 0.7 | PASS |
| 20Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 220 | 1.5 | PASS |
| 24Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 277 | 0.6 | PASS |
| 32Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 409 | 0.8 | PASS |
| 48Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 706 | 2.2 | PASS |
| 64Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1082 | 3.6 | PASS |
| 96Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1808 | 9.9 | PASS |
| 128Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 2690 | 23.4 | PASS |

- correctness/stability ceiling (all workloads): **128Q**
- runtime ceiling (median opt time sane): **128Q**
- memory ceiling (< 1024 MB peak): **128Q**

Verified-incorrect anywhere is a release blocker; 'policy-only' rows are stability-proven, not equivalence-proven (see README tier table).
