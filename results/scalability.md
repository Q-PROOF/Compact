# Scalability gauntlet (generated - do not edit)

Generated 2026-09-18T18:26:31.854532+00:00 | commit `53d00ca20e3f` | compactq 0.1.6 | referee limit 10q (independent Qiskit Operator); Clifford workloads are algebraically proven at any width.

| width | tests | crashes | timeouts | incorrect | policy viol | nondet | verified | median opt ms | peak MB | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 75 | 0.9 | PASS |
| 6Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 124 | 2.0 | PASS |
| 8Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 93 | 0.3 | PASS |
| 10Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 138 | 0.4 | PASS |
| 12Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 148 | 0.5 | PASS |
| 16Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 208 | 0.7 | PASS |
| 20Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 228 | 1.5 | PASS |
| 24Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 276 | 0.6 | PASS |
| 32Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 410 | 0.8 | PASS |
| 48Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 710 | 2.2 | PASS |
| 64Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1070 | 3.6 | PASS |
| 96Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1786 | 9.9 | PASS |
| 128Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 2660 | 23.4 | PASS |
| 256Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 28773 | 168.0 | PASS |

- correctness/stability ceiling (all workloads): **256Q**
- runtime ceiling (median ≤ 30 s, 'practical' class): **256Q**
- memory ceiling (< 1024 MB peak): **256Q**

Verified-incorrect anywhere is a release blocker; 'policy-only' rows are stability-proven, not equivalence-proven (see README tier table).
