# Independent reproduction harness (generated)

Generated 2026-09-19T04:51:21.444447+00:00 | commit `d9caa68dec81` | tools: compact, qiskit, pytket, cirq

| circuit | qubits | in 2q (u3+cx) | compact 2q | qiskit 2q | pytket 2q | cirq 2q | compact status | qiskit status | pytket status | cirq status |
|---|---:|---:| ---:|---:|---:|---:| ---|---|---|---|
| qft_4 | 4 | 18 | 18 | 12 | 12 | 18 | OK | INEQUIVALENT | INEQUIVALENT | OK |
| qft_6 | 6 | 39 | 39 | 30 | 30 | 39 | OK | INEQUIVALENT | INEQUIVALENT | OK |
| qft_8 | 8 | 68 | 68 | 56 | 56 | 68 | OK | INEQUIVALENT | INEQUIVALENT | OK |
| qaoa_4 | 4 | 12 | 12 | 12 | 12 | 12 | OK | OK | OK | OK |
| qaoa_6 | 6 | 20 | 20 | 20 | 20 | 20 | OK | OK | OK | OK |
| qaoa_8 | 8 | 28 | 28 | 28 | 28 | 28 | OK | OK | OK | OK |
| qpe_4 | 4 | 11 | 9 | 8 | 6 | ? | OK | INEQUIVALENT | INEQUIVALENT | error ValueError: operands could not be broadcast together with shapes (16,16) (8,8)  |
| qpe_6 | 6 | 23 | 14 | 20 | 11 | ? | OK | INEQUIVALENT | INEQUIVALENT | error ValueError: operands could not be broadcast together with shapes (64,64) (16,16)  |
| grover_3 | 3 | 24 | 24 | 24 | 14 | 32 | OK | OK | INEQUIVALENT | OK |
| grover_5 | 5 | 144 | 386 | 144 | 144 | 144 | OK | OK | OK | OK |
| heisenberg_4 | 4 | 36 | 18 | 36 | 18 | 18 | OK | OK | OK | OK |
| heisenberg_6 | 6 | 60 | 30 | 60 | 30 | 30 | OK | OK | OK | OK |
| clifford_5 | 5 | 20 | 19 | 19 | 17 | 19 | OK | OK | INEQUIVALENT | OK |
| cliffordt_5 | 5 | 0 | 0 | 0 | 0 | 0 | OK | OK | OK | OK |
| adder_4 | 4 | 14 | 14 | 14 | 14 | 16 | OK | OK | OK | OK |
| adder_6 | 6 | 26 | 26 | 26 | 26 | 30 | OK | OK | OK | OK |
| vqe_4 | 4 | 3 | 3 | 3 | 3 | 3 | OK | OK | OK | OK |
| vqe_6 | 6 | 5 | 5 | 5 | 5 | 5 | OK | OK | OK | OK |
| redundant_4 | 4 | 24 | 0 | 0 | 0 | ? | OK | OK | OK | error ValueError: operands could not be broadcast together with shapes (16,16) (8,8)  |
| redundant_6 | 6 | 40 | 0 | 0 | 0 | ? | OK | OK | OK | error ValueError: operands could not be broadcast together with shapes (64,64) (32,32)  |

## Summary

- **compact**: mean 2q cut 10.2%, median wall 21 ms
- **qiskit**: mean 2q cut 14.6%, median wall 4 ms
- **pytket**: mean 2q cut 25.0%, median wall 175 ms
- **cirq**: mean 2q cut 2.8%, median wall 51 ms
- compact vs qiskit 2q W/T/L: [2, 12, 1]
