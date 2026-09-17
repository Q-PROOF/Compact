"""Resource estimation: the Q#-ResourceEstimator-style accounting of
what a circuit costs at the fault-tolerance layer.

Exact counts over the circuit IR (no sampling, no heuristics):
  t_count    T + Tdg gates - the magic-state currency of FT hardware;
             each must be produced by a distillation factory
  t_depth    longest chain of T/Tdg layers (distillation latency bound)
  cx_count / two_qubit / gates / depth / num_qubits
  clifford_1q  everything else (virtually free in FT)

The FT cost model is deliberately NOT fabricated: lattice-surgery
distillation costs depend on code distance, factory layouts and
decoder hardware, so this module reports the exact gate-level
resource vector that FT cost models (e.g. direct lattice-surgery
compilation pipelines) consume as input.
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate

_T_GATE_NAMES = {"t", "tdg"}
_CLIFFORD_1Q = {"h", "x", "y", "z", "s", "sdg", "sx", "sxdg"}


def resource_estimate(circ: Circuit) -> dict:
    """Exact gate-level resource vector of `circ`."""
    t_count = 0
    cx = 0
    two_q = 0
    clifford_1q = 0
    non_clifford_1q = 0
    for g in circ.ops:
        if g.name in _T_GATE_NAMES:
            t_count += 1
        elif g.name == "cx":
            cx += 1
        elif len(g.qubits) == 2:
            two_q += 1
        elif g.name in _CLIFFORD_1Q:
            clifford_1q += 1
        elif g.name in ("rx", "ry", "rz", "p") and (not g.params
                or _is_clifford_angle(float(g.params[0]))):
            clifford_1q += 1
        else:
            non_clifford_1q += 1
    return {
        "num_qubits": circ.num_qubits,
        "gates": len(circ.ops),
        "two_qubit": circ.two_qubit_count(),
        "cx_count": cx,
        "other_2q": two_q,
        "depth": circ.depth(),
        "t_count": t_count,
        "clifford_1q": clifford_1q,
        "non_clifford_1q": non_clifford_1q,
    }


def t_depth(circ: Circuit) -> int:
    """Longest chain of T/Tdg gates along any wire (the distillation
    latency bound), computed as an ASAP T-layer count."""
    wire_layer: dict[int, int] = {}
    max_layer = 0
    for g in circ.ops:
        if g.name not in _T_GATE_NAMES:
            continue
        start = max((wire_layer.get(w, 0) for w in g.qubits), default=0)
        for w in g.qubits:
            wire_layer[w] = start + 1
        max_layer = max(max_layer, start + 1)
    return max_layer


def _is_clifford_angle(theta: float) -> bool:
    k = round(theta / (math.pi / 2))
    return abs(theta - k * math.pi / 2) < 1e-9


def rebase_cliffordt(circ: Circuit) -> Circuit:
    """Re-express Clifford-angle Rz/P rotations into the {H, S, T}
    Clifford+T frame exactly (up to global phase):
    rz(k*pi/4) = T^(k odd) . S^(k//2);  rz(k*pi/2) = S^k.
    Non-Clifford rotations pass through untouched (they are counted as
    non-Clifford by resource_estimate)."""
    out = []
    for g in circ.ops:
        if g.name in ("rz", "p") and g.params:
            theta = float(g.params[0])
            k = round(theta / (math.pi / 4))
            if abs(theta - k * (math.pi / 4)) < 1e-9:
                kk = k % 8
                for _ in range(kk // 2):
                    out.append(Gate("s", (), g.qubits))
                if kk % 2:
                    out.append(Gate("t", (), g.qubits))
                continue
        out.append(g)
    return Circuit(circ.num_qubits, out)
