"""Optional pytket interoperability (requires `pip install pytket`).

Loss-free circuit conversion between compactq's IR and pytket
``Circuit``s, mirroring :mod:`compactq.qiskit_bridge`:

    from compactq.pytket_bridge import from_pytket, to_pytket

Conventions handled for you:

* pytket rotations (Rx/Ry/Rz/CRz/TK1, ...) are measured in HALF-TURNS;
  compactq uses radians — converted exactly here (verified against
  tket's own matrices to 1e-13).
* ``cp`` is emitted as the exact Rz/CX controlled-phase decomposition
  (tket's ``CRz`` has a different phase convention and is expanded on
  import instead).
* Measurement, reset and classical operations raise
  UnsupportedCircuitError — compactq optimizes unitary circuits.
* Gates outside the conversion table are rebased with pytket's own
  ``AutoRebase`` onto the shared gate set; if that fails, the error is
  precise, never silent.
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate, ONE_Q_GATES, TWO_Q_GATES, MULTI_Q_GATES
from .errors import UnsupportedCircuitError

__all__ = ["to_pytket", "from_pytket"]

# pytket OpTypes expressible directly in compactq's IR
_1Q_TABLE = {
    "H": "h", "X": "x", "Y": "y", "Z": "z", "S": "s", "Sdg": "sdg",
    "T": "t", "Tdg": "tdg",
}
_HALF_TURN_1Q = {"Rx": "rx", "Ry": "ry", "Rz": "rz"}
_2Q_TABLE = {"CX": "cx", "CZ": "cz", "SWAP": "swap"}


def to_pytket(circ: Circuit):
    """Convert a compactq Circuit into a pytket Circuit (loss-free up to
    the global phase conventions both projects ignore)."""
    from pytket.circuit import Circuit as TCircuit

    tk = TCircuit(circ.num_qubits)
    for g in circ.ops:
        name, params, qs = g.name, list(g.params), list(g.qubits)
        if name in ("h", "x", "y", "z", "s", "sdg", "t", "tdg"):
            getattr(tk, {"h": "H", "x": "X", "y": "Y", "z": "Z", "s": "S",
                         "sdg": "Sdg", "t": "T", "tdg": "Tdg"}[name])(qs[0])
        elif name in ("rx", "ry", "rz"):
            getattr(tk, name.capitalize())(params[0] / math.pi, qs[0])
        elif name == "p":            # diag(1, e^{i th}) == Rz up to phase
            tk.Rz(params[0] / math.pi, qs[0])
        elif name == "sx":           # rx(pi/2) up to global phase
            tk.Rx(0.5, qs[0])
        elif name == "sxdg":
            tk.Rx(-0.5, qs[0])
        elif name in ("u", "u3"):    # rz(lam) . ry(theta) . rz(phi)
            lam, theta, phi = params[2], params[0], params[1]
            tk.Rz(lam / math.pi, qs[0])
            tk.Ry(theta / math.pi, qs[0])
            tk.Rz(phi / math.pi, qs[0])
        elif name == "cx":
            tk.CX(qs[0], qs[1])
        elif name == "cz":
            tk.CZ(qs[0], qs[1])
        elif name == "swap":
            tk.SWAP(qs[0], qs[1])
        elif name == "cp":           # exact: P(l/2)c P(l/2)t CX Rz(-l/2)t CX
            l = params[0]
            tk.Rz(l / (2 * math.pi), qs[0])
            tk.Rz(l / (2 * math.pi), qs[1])
            tk.CX(qs[0], qs[1])
            tk.Rz(-l / (2 * math.pi), qs[1])
            tk.CX(qs[0], qs[1])
        else:
            # ecr / iswap / mcx / mcp: expand into the core set first
            from .mcx import expand_gate
            ex = expand_gate(g)
            if ex is None:
                raise UnsupportedCircuitError(
                    f"no pytket conversion for gate {name!r}")
            tk.add_circuit(to_pytket(Circuit(circ.num_qubits, list(ex))))
    return tk


def from_pytket(tk) -> Circuit:
    """Convert a pytket Circuit into compactq's IR, loss-free.

    Raises UnsupportedCircuitError when the circuit contains measurement,
    reset, or classical operations, or gates neither of the two toolsets
    can rebase.  Barriers are skipped (no-ops)."""
    from pytket.passes import AutoRebase, DecomposeBoxes
    from pytket.circuit import OpType

    tk = tk.copy()
    DecomposeBoxes().apply(tk)
    _SUPPORTED = {OpType.H, OpType.X, OpType.Y, OpType.Z, OpType.S,
                  OpType.Sdg, OpType.T, OpType.Tdg, OpType.Rx, OpType.Ry,
                  OpType.Rz, OpType.CX, OpType.CZ, OpType.SWAP,
                  OpType.CCX, OpType.CSWAP, OpType.TK1, OpType.CRz,
                  OpType.Barrier}
    names = {c.op.type for c in tk}
    if not names <= _SUPPORTED:
        try:
            AutoRebase({OpType.H, OpType.X, OpType.S, OpType.T,
                        OpType.Rx, OpType.Ry, OpType.Rz, OpType.CX,
                        OpType.SWAP}).apply(tk)
        except Exception:
            bad = sorted(str(t) for t in names - _SUPPORTED)
            raise UnsupportedCircuitError(
                "pytket circuit uses gates compactq cannot rebase: "
                + ", ".join(bad))

    qmap = {q: i for i, q in enumerate(tk.qubits)}
    circ = Circuit(len(tk.qubits), [])
    from .mcx import expand_gate

    def half_turns(p):
        return float(p) * math.pi      # tket half-turns -> radians

    for cmd in tk:
        op = cmd.op
        t = op.type
        if t == OpType.Barrier:
            continue
        if t in (OpType.Measure, OpType.Reset, OpType.SetBits,
                 OpType.CopyBits, OpType.RangePredicate):
            raise UnsupportedCircuitError(
                "pytket circuit contains measurement/reset/classical "
                "operations; compactq optimizes unitary circuits")
        qs = tuple(qmap[q] for q in cmd.qubits)
        params = tuple(half_turns(p) for p in op.params)
        if t in _1Q_KEY:
            circ.append(Gate(_1Q_KEY[t], (), qs))
        elif t in _HT_KEY:
            circ.append(Gate(_HT_KEY[t], (params[0],), qs))
        elif t == OpType.CX:
            circ.append(Gate("cx", (), qs))
        elif t == OpType.CZ:
            circ.append(Gate("cz", (), qs))
        elif t == OpType.SWAP:
            circ.append(Gate("swap", (), qs))
        elif t == OpType.TK1:        # TK1(a,b,c) = Rz(a).Rx(b).Rz(c)
            a, b, c = params
            circ.append(Gate("rz", (c,), qs))
            circ.append(Gate("rx", (b,), qs))
            circ.append(Gate("rz", (a,), qs))
        elif t == OpType.CCX:        # toffoli -> native mcx (expanded)
            circ.append(Gate("mcx", (), qs))
        elif t == OpType.CSWAP:      # cx(t,c2) . ccx . cx(t,c2)
            c1, c2, tgt = qs
            circ.append(Gate("cx", (), (tgt, c2)))
            circ.append(Gate("mcx", (), (c1, c2, tgt)))
            circ.append(Gate("cx", (), (tgt, c2)))
        elif t == OpType.CRz:
            # tket CRz(p) = P(-pi p / 2) on the control  x  CP(+pi p):
            # phase(-pi p/2) when control=1, CP adds +pi p on |control=1,
            # target=1> — both diagonal, so order is irrelevant
            circ.append(Gate("p", (-0.5 * params[0],), (qs[0],)))
            circ.append(Gate("cp", (params[0],), qs))
        else:                        # pragma: no cover (rebased above)
            g = Gate(str(t), (), qs) if not op.params else \
                Gate(str(t), params, qs)
            ex = expand_gate(g)
            if ex is None:
                raise UnsupportedCircuitError(
                    f"pytket gate {t} not convertible")
            circ.ops.extend(ex)
    # multi-qubit gates on repeated wires are impossible from tket, but
    # mcx/mcp re-validate through Gate construction anyway
    return circ


_1Q_KEY = {}
_HT_KEY = {}
try:  # populate the tables without hard-failing on pytket version drift
    from pytket.circuit import OpType as _OT
    for _ot, _nm in _1Q_TABLE.items():
        _1Q_KEY[getattr(_OT, _ot)] = _nm
    for _ot, _nm in _HALF_TURN_1Q.items():
        _HT_KEY[getattr(_OT, _ot)] = _nm
except Exception:  # pragma: no cover - pytket not installed
    pass
