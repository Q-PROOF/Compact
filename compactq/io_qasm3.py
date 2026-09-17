"""OpenQASM 3 export for the compactq gate set."""
from __future__ import annotations

from .circuit import Circuit


def to_qasm3(circ: Circuit, name: str = "q") -> str:
    lines = ["OPENQASM 3.0;", 'include "stdgates.inc";', f"qubit[{circ.num_qubits}] {name};"]
    for g in circ.ops:
        qs = ", ".join(f"{name}[{q}]" for q in g.qubits)
        if g.name == "p":
            lines.append(f"p({g.params[0]:.17g}) {qs};")
        elif g.params:
            lines.append(f"{g.name}({g.params[0]:.17g}) {qs};")
        else:
            lines.append(f"{g.name} {qs};")
    return "\n".join(lines) + "\n"
