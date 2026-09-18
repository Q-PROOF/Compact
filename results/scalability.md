# Scalability gauntlet (generated - do not edit)

Generated 2026-09-18T17:08:48.024466+00:00 | commit `fba3e40b9d22` | compactq 0.1.3 | referee limit 10q (independent Qiskit Operator); Clifford workloads are algebraically proven at any width.

| width | tests | crashes | timeouts | incorrect | policy viol | nondet | verified | median opt ms | peak MB | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 76 | 0.9 | PASS |
| 6Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 124 | 2.0 | PASS |
| 8Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 92 | 0.3 | PASS |
| 10Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 140 | 0.4 | PASS |
| 12Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 146 | 0.5 | PASS |
| 16Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 210 | 0.7 | PASS |
| 20Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 224 | 1.5 | PASS |
| 24Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 294 | 0.6 | PASS |
| 32Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 442 | 0.8 | PASS |
| 48Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 721 | 2.2 | PASS |
| 64Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1094 | 3.6 | PASS |
| 96Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1816 | 9.9 | PASS |
| 128Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 2716 | 23.4 | PASS |
| 256Q | 8 | 0 | 1 | 0 | 0 | 0 | 2 | 34252 | 168.0 | TEST |

- correctness/stability ceiling (all workloads): **128Q**
- runtime ceiling (median opt time sane): **128Q**
- memory ceiling (< 1024 MB peak): **128Q**

Verified-incorrect anywhere is a release blocker; 'policy-only' rows are stability-proven, not equivalence-proven (see README tier table).
