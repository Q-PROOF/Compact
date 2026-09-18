# Scalability gauntlet (generated - do not edit)

Generated 2026-09-18T17:45:24.534796+00:00 | commit `fe4288d32e04` | compactq 0.1.4 | referee limit 10q (independent Qiskit Operator); Clifford workloads are algebraically proven at any width.

| width | tests | crashes | timeouts | incorrect | policy viol | nondet | verified | median opt ms | peak MB | status |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 78 | 0.9 | PASS |
| 6Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 123 | 2.0 | PASS |
| 8Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 100 | 0.3 | PASS |
| 10Q | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 137 | 0.4 | PASS |
| 12Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 142 | 0.5 | PASS |
| 16Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 208 | 0.7 | PASS |
| 20Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 232 | 1.5 | PASS |
| 24Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 250 | 0.6 | PASS |
| 32Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 434 | 0.8 | PASS |
| 48Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 735 | 2.2 | PASS |
| 64Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1136 | 3.6 | PASS |
| 96Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 1720 | 9.9 | PASS |
| 128Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 2694 | 23.4 | PASS |
| 256Q | 8 | 0 | 0 | 0 | 0 | 0 | 2 | 33968 | 168.0 | PASS |

- correctness/stability ceiling (all workloads): **256Q**
- runtime ceiling (median opt time sane): **256Q**
- memory ceiling (< 1024 MB peak): **256Q**

Verified-incorrect anywhere is a release blocker; 'policy-only' rows are stability-proven, not equivalence-proven (see README tier table).
