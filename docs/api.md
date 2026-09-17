# API

```python
import compactq
compactq.optimize(circ, verify=True)          # safe pipeline, proven exact
compactq.optimize_deep(circ)                  # KAK cascade for dense circuits
compactq.optimize_search(circ)                # multi-pipeline search, best kept
compactq.approximate(circ, min_fidelity=0.99) # bounded-fidelity 2q reduction
```

## Hardware-aware

```python
from compactq.target import Target, optimize_for
t = Target(cx_fidelity={(0, 1): 0.999, (1, 0): 0.98},
           default_cx_fidelity=0.99,
           single_qubit_fidelity=0.9999)
best, est_fidelity = optimize_for(circ, t)
```

## Clifford

```python
compactq.is_clifford(circ)          # bool
compactq.clifford_equal(a, b)       # exact equality proof, any qubit count
```

## IO

```python
compactq.from_qasm(text) / compactq.to_qasm(circ)   # OpenQASM 2.0
compactq.from_qasm3(text)                       # OpenQASM 3 common subset
compactq.to_qasm3(circ)
```

`from_qasm` understands the extended qelib1 set (swap, cswap, sx, sxdg, cy,
crz, cu1/cp, rzz, rxx, ccx, u/u3/u2/u1), user gate definitions, multiple
registers, and register-wide operands.  `mcx`/`mcp`/`c3x`/`c4x` import as
native multi-controlled gates.

## Hardware-aware passes

```python
from compactq.hardware import route, flip_cx, translate_1q_to_rz_sx_x
routed, final_map = route(circ, coupling=[{0, 1}, {1, 2}])
```

## Optional Rust kernels

```bash
cd native && pip install maturin && maturin build --release -o dist
pip install dist/compactq_native-*.whl
```

Raises the exact-verification ceiling to 8 qubits and accelerates the KAK
hot path (~4x).  compactq auto-detects the wheel and silently falls back.
