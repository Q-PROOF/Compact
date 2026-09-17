"""Optional Cirq interoperability (export; requires `pip install cirq`)."""
from __future__ import annotations

import math

from .circuit import Circuit


def to_cirq(circ: Circuit):
    import cirq

    qubits = [cirq.LineQubit(i) for i in range(circ.num_qubits)]
    ops = []
    for g in circ.ops:
        qs = [qubits[q] for q in g.qubits]
        n = g.name
        simple = {"h": cirq.H, "x": cirq.X, "y": cirq.Y, "z": cirq.Z,
                  "s": cirq.S, "t": cirq.T, "sx": cirq.SQRT_X}
        if n in simple:
            ops.append(simple[n](*qs))
        elif n == "sdg":
            ops.append(cirq.S(*qs) ** -1)
        elif n == "tdg":
            ops.append(cirq.T(*qs) ** -1)
        elif n in ("rx", "ry", "rz"):
            ops.append(getattr(cirq, n)(math.degrees(g.params[0]))(*qs))
        elif n == "p":
            ops.append(cirq.ZPowGate(exponent=g.params[0] / math.pi)(*qs))
        elif n == "cx":
            ops.append(cirq.CNOT(*qs))
        elif n == "cz":
            ops.append(cirq.CZ(*qs))
        elif n == "swap":
            ops.append(cirq.SWAP(*qs))
        else:
            raise ValueError(f"cirq bridge: unsupported gate {n!r}")
    return cirq.Circuit(ops)
