# Feynman-corpus verified optimization (generated)

Generated 2026-09-24T04:22:07.091087+00:00 | corpus: meamy/feynman benchmarks/qasm (42 parsed, 2 skipped)

compactq optimize_search (proof included) — 18/42 verified in-product

| circuit | q | in gates | out gates | in 2q | out 2q | wall ms | proof | referee |
|---|---:|---:|---:|---:|---:|---:|---|---|
| adder_8 | 24 | 1128 | 743 | 409 | 340 | 7813.5 | prover_unavailable | in-product |
| barenco_tof_10 | 19 | 578 | 340 | 192 | 160 | 471.8 | prover_unavailable | in-product |
| barenco_tof_3 | 5 | 76 | 42 | 24 | 20 | 50.4 | full_unitary | None |
| barenco_tof_4 | 7 | 146 | 82 | 48 | 40 | 115.1 | full_unitary | None |
| barenco_tof_5 | 9 | 218 | 122 | 72 | 60 | 155.7 | dd_full_unitary | in-product |
| csla_mux_3 | 15 | 210 | 137 | 80 | 63 | 101.7 | randomized_sampling | in-product |
| csum_mux_9 | 30 | 532 | 313 | 168 | 140 | 678.7 | prover_unavailable | in-product |
| cycle_17_3 | - | - | - | - | - | - | skipped (ValueError: 'cx' acts on a repeated qubit: (28, 28)) | - |
| gf2^10_mult | 30 | 1747 | 1083 | 609 | 509 | 479.6 | prover_unavailable | in-product |
| gf2^128_mult | 384 | 279419 | 213881 | 98685 | 98685 | 6366.4 | prover_unavailable | in-product |
| gf2^16_mult | 48 | 4459 | 3433 | 1581 | 1581 | 103.4 | prover_unavailable | in-product |
| gf2^256_mult | 768 | 1115899 | 853753 | 393981 | 393981 | 24253.6 | prover_unavailable | in-product |
| gf2^32_mult | 96 | 17658 | 13560 | 6268 | 6268 | 403.1 | prover_unavailable | in-product |
| gf2^4_mult | 12 | 289 | 163 | 99 | 83 | 152.6 | randomized_sampling | in-product |
| gf2^5_mult | 15 | 447 | 249 | 154 | 129 | 249.0 | randomized_sampling | in-product |
| gf2^64_mult | 192 | 70075 | 53689 | 24765 | 24765 | 1625.0 | prover_unavailable | in-product |
| gf2^6_mult | 18 | 639 | 353 | 221 | 185 | 448.5 | prover_unavailable | in-product |
| gf2^7_mult | 21 | 865 | 479 | 300 | 251 | 594.6 | prover_unavailable | in-product |
| gf2^8_mult | 24 | 1139 | 641 | 405 | 341 | 836.4 | prover_unavailable | in-product |
| gf2^9_mult | 27 | 1419 | 869 | 494 | 413 | 323.3 | prover_unavailable | in-product |
| grover_5 | 9 | 1023 | 521 | 288 | 240 | 621.8 | randomized_sampling | in-product |
| ham15-high | 20 | 6712 | 4917 | 2149 | 2149 | 421.9 | prover_unavailable | in-product |
| ham15-low | 17 | 535 | 349 | 236 | 212 | 304.1 | dd_full_unitary | in-product |
| ham15-med | 17 | 1600 | 1033 | 534 | 452 | 364.8 | prover_unavailable | in-product |
| hwb10 | 16 | 91642 | 70063 | 35170 | 34814 | 7198.4 | prover_unavailable | in-product |
| hwb11 | 15 | 256181 | 195549 | 98023 | 97621 | 20219.7 | prover_unavailable | in-product |
| hwb12 | 20 | 514412 | 389740 | 191803 | 191411 | 43875.2 | prover_unavailable | in-product |
| hwb6 | 7 | 319 | 196 | 116 | 100 | 216.1 | full_unitary | None |
| hwb8 | 12 | 18220 | 14148 | 7129 | 7099 | 1056.4 | prover_unavailable | in-product |
| mod5_4 | 5 | 79 | 46 | 28 | 24 | 35.3 | full_unitary | None |
| mod_adder_1024 | 28 | 5425 | 3957 | 1720 | 1702 | 423.3 | prover_unavailable | in-product |
| mod_adder_1048576 | - | - | - | - | - | - | skipped (ValueError: 'cx' acts on a repeated qubit: (48, 48)) | - |
| mod_mult_55 | 9 | 147 | 93 | 48 | 42 | 131.5 | dd_full_unitary | in-product |
| mod_red_21 | 11 | 346 | 191 | 105 | 88 | 190.6 | dd_full_unitary | in-product |
| qcla_adder_10 | 36 | 657 | 410 | 233 | 189 | 618.3 | prover_unavailable | in-product |
| qcla_com_7 | 24 | 559 | 332 | 186 | 151 | 411.5 | prover_unavailable | in-product |
| qcla_mod_7 | 26 | 1120 | 671 | 382 | 315 | 850.7 | prover_unavailable | in-product |
| qft_4 | 5 | 187 | 86 | 46 | 39 | 88.1 | full_unitary | None |
| rc_adder_6 | 14 | 244 | 157 | 93 | 76 | 112.3 | dd_full_unitary | in-product |
| tof_10 | 19 | 323 | 189 | 102 | 85 | 313.0 | prover_unavailable | in-product |
| tof_3 | 5 | 57 | 33 | 18 | 15 | 35.9 | full_unitary | None |
| tof_4 | 7 | 95 | 55 | 30 | 25 | 71.2 | full_unitary | None |
| tof_5 | 9 | 133 | 79 | 42 | 35 | 81.9 | dd_full_unitary | in-product |
| vbe_adder_3 | 10 | 190 | 114 | 70 | 56 | 77.9 | dd_full_unitary | in-product |
