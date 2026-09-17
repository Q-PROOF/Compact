# CLI

```bash
python -m compactq INPUT.qasm -o OUTPUT.qasm --stats
compactq INPUT.qasm --stats              # installed console script
```

Options:

| flag | meaning |
|---|---|
| `-o/--output` | write result to a file (default stdout) |
| `--stats` | print before/after statistics to stderr |
| `--approx FIDELITY` | approximate mode with per-block fidelity floor |
| `--no-verify` | skip the exact whole-circuit proof |

Benchmark suite:

```bash
python -m compactq.bench            # full comparison vs Qiskit (if installed)
python -m compactq.bench --quick    # 3-circuit CI guard
compactq-bench --out results.md
```
