# Benchmarks

All numbers generated on this repo (Python 3.11, qiskit 2.5.2). Fairness:
same input QASM, same basis set, fixed seeds. `gates / 2q / depth`.
Circuits with classical control flow or mid-circuit measurement are skipped
with a stated reason (a unitary optimizer cannot represent them).

See the repository README for the current full tables (synthetic suite +
QASMBench unitary cores vs Qiskit L3 and pytket).

Highlights vs Qiskit L3 (full tables + competitor landscape in the README):

| circuit | compactq | qiskit L3 |
|---|---|---|
| error_correctiond3_n5 | **23 / 8 / 13** | 93 / 35 / 66 |
| pea_n5 | **34 / 10 / 21** | 51 / 17 / 32 |
| qaoa_n6 | **114 / 36 / 47** | 171 / 36 / 64 |
| ising_n10 | **166 / 49 / 29** | 260 / 90 / 46 |
| hhl_n7 | **191 / 72 / 128** | 253 / 92 / 169 |
| qft_n4 | **20 / 6 / 10** | 34 / 12 / 20 |
| adder_n10 | **110 / 57 / 95** | 137 / 65 / 99 |
| adder_n4 | **16 / 7 / 9** | 23 / 10 / 11 |

Reproduce with `python scripts/realbench.py` (qiskit + pytket optional).
