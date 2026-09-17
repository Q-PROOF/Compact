# Benchmark results (generated - do not edit)

Generated 2026-09-17T18:26:29.128531+00:00 | commit `9dc08b12d11c` | compactq 0.1.2 | qiskit 2.5.2 | native yes 0.16.4

| circuit | q | Compact T/2q/depth | Qiskit L3 T/2q/depth | proof | referee fid |
|---|---|---|---|---|---|
| adder_n10 | 10 | 110/57/95 | 137/65/99 | unverified | - |
| adder_n4 | 4 | 16/7/9 | 23/10/11 | exact-unitary | 1.0 |
| basis_change_n3 | 3 | 34/10/22 | 44/10/27 | exact-unitary | 1.0 |
| basis_test_n4 | 4 | 34/12/15 | 34/6/12 | exact-unitary | 1.0 |
| basis_trotter_n4 | 4 | 773/240/352 | 798/179/358 | exact-unitary | 1.0 |
| bb84_n8 | - | skipped: non-unitary (mid-circuit measurement) | - | - | - |
| bell_n4 | 4 | 18/5/7 | 29/5/12 | exact-unitary | 1.0 |
| cat_state_n4 | 4 | 4/3/4 | 4/3/4 | exact-unitary | 1.0 |
| deutsch_n2 | 2 | 4/1/3 | 4/1/3 | exact-unitary | 1.0 |
| dnn_n2 | 2 | 12/3/8 | 20/3/13 | exact-unitary | 1.0 |
| dnn_n8 | 8 | 216/64/37 | 340/64/60 | exact-unitary | 1.0 |
| error_correctiond3_n5 | 5 | 23/8/13 | 91/35/65 | exact-unitary | 1.0 |
| fredkin_n3 | 3 | 19/8/11 | 19/8/11 | exact-unitary | 1.0 |
| grover_n2 | 2 | 7/2/5 | 7/2/5 | exact-unitary | 1.0 |
| hhl_n10 | - | skipped: parse: QASM2ParseError: "<input>:186801,8: 'q' is not defined in this scope" | - | - | - |
| hhl_n7 | 7 | 191/72/128 | 255/92/168 | exact-unitary | 1.0 |
| hs4_n4 | 4 | 12/4/5 | 12/4/5 | exact-unitary | 1.0 |
| inverseqft_n4 | - | skipped: classical control flow | - | - | - |
| ipea_n2 | - | skipped: classical control flow | - | - | - |
| ising_n10 | 10 | 166/49/29 | 260/90/46 | unverified | - |
| iswap_n2 | 2 | 7/2/5 | 8/2/6 | exact-unitary | 1.0 |
| linearsolver_n3 | 3 | 11/4/9 | 15/4/10 | exact-unitary | 1.0 |
| lpn_n5 | 5 | 7/2/4 | 7/2/4 | exact-unitary | 1.0 |
| pea_n5 | 5 | 34/10/21 | 49/17/29 | exact-unitary | 1.0 |
| qaoa_n3 | - | skipped: non-unitary (mid-circuit measurement) | - | - | - |
| qaoa_n6 | 6 | 114/36/47 | 165/36/63 | exact-unitary | 1.0 |
| qec_en_n5 | 5 | 23/10/15 | 23/10/15 | exact-unitary | 1.0 |
| qec_sm_n5 | - | skipped: classical control flow | - | - | - |
| qft_n4 | 4 | 20/6/10 | 34/12/20 | exact-unitary | 1.0 |
| qpe_n9 | - | skipped: non-unitary (mid-circuit measurement) | - | - | - |
| qrng_n4 | 4 | 4/0/1 | 4/0/1 | exact-unitary | 1.0 |
| quantumwalks_n2 | 2 | 8/2/5 | 20/3/13 | exact-unitary | 0.999999999972 |
| sat_n7 | 7 | 125/52/71 | 158/60/86 | exact-unitary | 1.0 |
| shor_n5 | - | skipped: classical control flow | - | - | - |
| simon_n6 | 6 | 30/12/22 | 43/14/27 | exact-unitary | 1.0 |
| teleportation_n3 | 3 | 5/2/4 | 6/2/4 | exact-unitary | 1.0 |
| toffoli_n3 | 3 | 14/5/10 | 18/6/12 | exact-unitary | 1.0 |
| variational_n4 | 4 | 29/8/14 | 45/8/19 | exact-unitary | 1.0 |
| vqe_n4 | 4 | 25/9/11 | 46/9/18 | exact-unitary | 1.0 |
| vqe_uccsd_n4 | - | skipped: parse: QASM2ParseError: "<input>:225,8: 'q' is not defined in this scope" | - | - | - |
| vqe_uccsd_n6 | - | skipped: parse: QASM2ParseError: "<input>:2286,8: 'q' is not defined in this scope" | - | - | - |
| vqe_uccsd_n8 | - | skipped: parse: QASM2ParseError: "<input>:10813,8: 'q' is not defined in this scope" | - | - | - |
| wstate_n3 | 3 | 18/6/12 | 20/6/13 | exact-unitary | 1.0 |
