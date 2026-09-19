# Full-transpilation comparison (generated)

Generated 2026-09-19T08:12:23.742532+00:00 | commit `371ce35e9c78` | trivial layout both arms, sabre routing (qiskit), swaps included, both outputs refereed

| topology | circuit | compact 2q | qiskit 2q | compact depth | qiskit depth | fid c/q | winner |
|---|---|---:|---:|---:|---:|---|---|
| line_5 | ghz_5 | 4 | 4 | 19 | 14 | 1.000000/1.000000 | tie |
| line_5 | qaoa_ring_5 | 16 | 16 | 98 | 50 | 1.000000/1.000000 | tie |
| line_5 | qft_5 | 80 | 36 | 182 | 110 | 1.000000/0.250000 | qiskit |
| line_5 | heisenberg_5 | 24 | 24 | 172 | 69 | 1.000000/1.000000 | tie |
| line_5 | vqe_5 | 12 | 12 | 50 | 36 | 1.000000/1.000000 | tie |
| grid_4x4 | ghz_5 | 22 | 11 | 52 | 21 | 1.000000/0.125000 | qiskit |
| grid_4x4 | qaoa_ring_5 | 46 | 25 | 171 | 66 | 1.000000/0.125000 | qiskit |
| grid_4x4 | qft_5 | 92 | 38 | 243 | 102 | 1.000000/0.250000 | qiskit |
| grid_4x4 | heisenberg_5 | 66 | 43 | 283 | 119 | 1.000000/0.125000 | qiskit |
| grid_4x4 | vqe_5 | 54 | 31 | 131 | 69 | 1.000000/0.125000 | qiskit |
| grid_4x4 | ghz_16 | 87 | 44 | 119 | 83 | unverified-referee | qiskit |
| grid_4x4 | qaoa_8_on_grid | 58 | 40 | 202 | 91 | 1.000000/0.250000 | qiskit |
| grid_4x4 | qft_8_on_grid | 164 | 92 | 426 | 209 | 1.000000/0.015625 | qiskit |

post-routing 2q: compact wins 0 / ties 4 / losses 9 (losses reported, never hidden)

