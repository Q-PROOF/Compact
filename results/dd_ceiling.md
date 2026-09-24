# DD prover ceiling — structured families (generated)

Generated 2026-09-23T15:21:23.205925+00:00 | commit `ccde0ec41dd8` | budget 400000 nodes | normalization: max-weight (QMDD standard, v0.2.4)

| family | width | nodes | wall ms | status |
|---|---:|---:|---:|---|
| qft | 8 | 255 | 61 | proven |
| qft | 10 | 1023 | 244 | proven |
| qft | 12 | 4095 | 1083 | proven |
| qft | 14 | 16383 | 4756 | proven |
| qft | 16 | 65535 | 22683 | proven |
| qft | 18 | 262143 | 110703 | proven |
| qft | 20 | — | — | overflow (loud decline) |
| trotter_ring | 8 | 26 | 64 | proven |
| trotter_ring | 10 | 34 | 51 | proven |
| trotter_ring | 12 | 42 | 96 | proven |
| trotter_ring | 14 | 50 | 90 | proven |
| trotter_ring | 16 | 58 | 106 | proven |
| trotter_ring | 18 | 66 | 233 | proven |
| trotter_ring | 20 | 74 | 240 | proven |
| trotter_ring | 22 | 82 | 284 | proven |
| trotter_ring | 24 | 90 | 324 | proven |
| trotter_ring | 26 | 98 | 342 | proven |
| trotter_ring | 28 | 106 | 434 | proven |
| trotter_ring | 30 | 114 | 469 | proven |
| ghz | 8 | 15 | 4 | proven |
| ghz | 10 | 19 | 4 | proven |
| ghz | 12 | 23 | 5 | proven |
| ghz | 14 | 27 | 42 | proven |
| ghz | 16 | 31 | 26 | proven |
| ghz | 18 | 35 | 28 | proven |
| ghz | 20 | 39 | 16 | proven |
| ghz | 22 | 43 | 19 | proven |
| ghz | 24 | 47 | 25 | proven |
| ghz | 26 | 51 | 31 | proven |
| ghz | 28 | 55 | 42 | proven |
| ghz | 30 | 59 | 35 | proven |
| grover_mcx | 8 | 15 | 4732 | proven |
| grover_mcx | 10 | 19 | 39310 | proven |
| grover_mcx | 12 | 23 | 284075 | proven |
| grover_mcx | 14 | — | — | overflow (loud decline) |

## v0.2.3 -> v0.2.4 deltas

- **Correctness**: mcx chains with >= 7 controls produced a false INEQUIVALENT self-comparison verdict under v0.2.3 (self-fidelity 0.992 at 8q); fixed by max-weight normalization, with a self-consistency norm guard that converts any residual corruption into a loud decline.
- **Node counts**: unchanged for QFT (2^n - 1 — normalization does not change sharing); the honest node-budget ceiling for QFT-family proofs stays at 18q proven / 20q declined at the 400k default. A Rust DD kernel remains the roadmap item for 20-25q on this family; CX+diagonal circuits already prove algebraically at any width.
- **Grover-MCX family**: proven ceiling 10q -> 12q (v0.2.3 overflowed at 12q); 14q declines loudly.
- Trotter/GHZ families prove at every tested width to 30q, unchanged.
