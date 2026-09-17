# compactq

**Zero-dependency quantum circuit optimizer.** Smaller, shallower circuits —
exact by construction, verified before return, and fast.

```python
import compactq
from compactq.benchmarks import qft

opt = compactq.optimize_search(qft(4))
print(opt.stats())
```

`pip install compactq` — no dependencies, no account, no cloud. Pure Python >= 3.9.

## Why

Every gate you remove from a quantum circuit removes noise. compactq takes a
circuit and returns an **equivalent** one that is smaller (fewer gates),
shallower (lower depth), with priority on cutting 2-qubit gates — the
dominant error source on today's hardware.

## What's inside

- **Peephole folding** — maximal 1q runs resynthesized to <= 3 canonical
  gates, with a final single-`u3` fold for gate-count polish.
- **Commutation engine** — CX/CZ cancellation across provably-commuting
  gates, sliding passes, SWAP templates, generalized diagonal sliding.
- **CP engine** — `CX·RZ_t(θ)·CX = RZ_c(θ)·RZ_t(θ)·CP(−2θ)`: halves the CX
  count of QAOA/QFT/PEA-style phase ladders.
- **Pure-Python KAK/Weyl synthesis** — minimal-CX two-qubit synthesis,
  oracle-validated against Qiskit at ~1e-15.
- **Clifford resynthesis** — tableau-proven AG synthesis; Clifford blocks
  are proven exact at *any* qubit count.
- **Phase-polynomial pass** — parity-network re-synthesis of diagonal cores.
- **Approximate mode** — trade a bounded, measured amount of fidelity for
  fewer 2-qubit gates.
- **Hardware-aware objective** — optimize for estimated infidelity on *your*
  machine's error rates, not abstract gate counts.
- **Optional Rust kernels** — verified acceleration; auto-detected, silent
  pure-Python fallback.

## Quick links

- [Correctness model](correctness.md)
- [Benchmarks](benchmarks.md)
- [API reference](api.md)
- [CLI usage](cli.md)
