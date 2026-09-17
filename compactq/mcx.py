"""Multi-controlled gates: exact expansion into the core gate set.

MCP(theta) with controls c1..ck and target t is the diagonal gate
exp(i theta |1..1><1..1|) over the (k+1)-wire register.  Using the parity
expansion of the all-ones projector,

    |1..1><1..1| = (1/2^m) sum_{S subset W} (-1)^{|S|} Z_S,   m = k+1,

the gate decomposes exactly into RZ parity rotations on every non-empty
wire subset S - each realized by a CX ladder + RZ + CX ladder via the
parity network in compactq.parity.  MCX = H(t) . MCP(pi) . H(t).

Every expansion is unit-tested against qiskit's Operator for 1-3 controls.
"""
from __future__ import annotations

import math

from .circuit import Circuit, Gate


def expand_mcp(controls: tuple, target: int, theta: float) -> list[Gate]:
    """Exact cx+rz expansion of MCP(theta) (up to global phase)."""
    wires = list(controls) + [target]
    m = len(wires)
    out: list[Gate] = []
    for mask in range(1, 1 << m):
        support = 0
        for w in range(m):
            if (mask >> w) & 1:
                support |= 1 << wires[w]
        signs = bin(mask).count("1")
        rz_angle = theta * ((-1) ** (signs + 1)) / (2 ** (m - 1))
        if abs(rz_angle) < 1e-13:
            continue
        # realize exp(-i rz_angle Z_support / 2) with a CX ladder:
        # pick the lowest wire of the support as the phase wire, CNOT the
        # other support wires into it, rz, then undo the ladder.
        bits = [wires[w] for w in range(m) if (mask >> w) & 1]
        phase_wire = bits[0]
        ladder = [(b, phase_wire) for b in bits[1:]]
        for c, t in ladder:
            out.append(Gate("cx", (), (c, t)))
        out.append(Gate("rz", (rz_angle,), (phase_wire,)))
        for c, t in reversed(ladder):
            out.append(Gate("cx", (), (c, t)))
    return out


def expand_mcx(controls: tuple, target: int) -> list[Gate]:
    """Exact expansion of MCX: H(target) . MCP(pi) . H(target)."""
    out = [Gate("h", (), (target,))]
    out.extend(expand_mcp(controls, target, math.pi))
    out.append(Gate("h", (), (target,)))
    return out


def expand_gate(g: Gate) -> list[Gate] | None:
    """Expand an mcx/mcp Gate into core gates, or None if not one."""
    if g.name == "mcx":
        *ctrls, t = g.qubits
        if len(ctrls) < 1:
            return [Gate("x", (), (t,))]
        return expand_mcx(tuple(ctrls), t)
    if g.name == "mcp":
        *ctrls, t = g.qubits
        if len(ctrls) < 1:
            return [Gate("p", (g.params[0],), (t,))]
        return expand_mcp(tuple(ctrls), t, g.params[0])
    return None


def expand_circuit(circ: Circuit) -> Circuit:
    """Return a copy of the circuit with all mcx/mcp gates expanded."""
    if not any(g.name in ("mcx", "mcp") for g in circ.ops):
        return circ
    out = []
    for g in circ.ops:
        exp = expand_gate(g)
        if exp is None:
            out.append(g)
        else:
            out.extend(exp)
    return Circuit(circ.num_qubits, out)
