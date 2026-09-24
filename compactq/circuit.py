"""Minimal, zero-dependency quantum circuit intermediate representation (IR).

A Circuit is an ordered list of Gates acting on ``num_qubits`` qubits.
Gate semantics follow the OpenQASM 2 / qelib1 conventions (little-endian
qubit indexing, RZ defined as diag(e^{-i t/2}, e^{+i t/2})).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Sequence

GateName = str
Params = tuple
Qubits = tuple

ONE_Q_GATES = frozenset({"h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p", "sx", "sxdg", "u", "u3"})
TWO_Q_GATES = frozenset({"cx", "cz", "swap", "cp", "ecr", "iswap"})
MULTI_Q_GATES = frozenset({"mcx", "mcp"})  # k controls + target; expanded on use


@dataclass(frozen=True)
class Gate:
    """A single quantum gate application."""

    name: GateName
    params: Params
    qubits: Qubits

    def __post_init__(self) -> None:
        n = len(self.qubits)
        if n >= 2 and len(set(self.qubits)) != n:
            # e.g. cx(0,0) or mcx with control == target: no well-defined
            # unitary, and silently keeping it would make every later
            # proof meaningless — reject at construction, like unknown names
            raise ValueError(
                f"{self.name!r} acts on a repeated qubit: {self.qubits}")
        if n == 1 and self.name not in ONE_Q_GATES:
            raise ValueError(f"unknown 1-qubit gate: {self.name!r}")
        if n >= 2 and self.name in MULTI_Q_GATES:
            if n < 2:
                raise ValueError(f"{self.name} needs at least one control")
            if self.name == "mcp" and len(self.params) != 1:
                raise ValueError("mcp requires exactly one parameter")
            if self.name == "mcx" and self.params:
                raise ValueError("mcx takes no parameters")
            return
        if n == 2 and self.name not in TWO_Q_GATES:
            raise ValueError(f"unknown 2-qubit gate: {self.name!r}")
        if len(self.params) and n != 1 and self.name != "cp":
            raise ValueError("parameterised gates must act on exactly one qubit")
        if self.name == "cp" and len(self.params) != 1:
            raise ValueError("cp requires exactly one parameter")


@dataclass
class Circuit:
    """An ordered sequence of gates on ``num_qubits`` qubits."""

    num_qubits: int
    ops: list

    # ------------------------------------------------------------------ core
    def copy(self) -> "Circuit":
        return Circuit(self.num_qubits, list(self.ops))

    def append(self, gate: Gate) -> None:
        for q in gate.qubits:
            if not 0 <= q < self.num_qubits:
                raise ValueError(f"qubit index {q} out of range for {self.num_qubits}-qubit circuit")
        self.ops.append(gate)

    def __len__(self) -> int:
        return len(self.ops)

    def __iter__(self):
        return iter(self.ops)

    # ------------------------------------------------------------- metrics
    def gate_counts(self) -> Counter:
        return Counter(g.name for g in self.ops)

    def two_qubit_count(self) -> int:
        return sum(1 for g in self.ops if len(g.qubits) == 2)

    def cx_equivalent_count(self) -> int:
        """Level-B hardware-cost metric: 2-qubit gates counted in
        CX-equivalents (a SWAP decomposes into 3 CX; every other 2-qubit
        gate counts 1).  The logical count is `two_qubit_count`."""
        return sum(3 if g.name == "swap" else 1
                   for g in self.ops if len(g.qubits) == 2)

    def depth(self) -> int:
        """Standard gate-depth (ASAP layering)."""
        if not self.ops:
            return 0
        frontier = [0] * self.num_qubits
        for g in self.ops:
            t = max(frontier[q] for q in g.qubits) + 1
            for q in g.qubits:
                frontier[q] = t
        return max(frontier)

    # -------------------------------------------------------------- display
    def stats(self) -> str:
        c = self.gate_counts()
        parts = [f"{v} {k}" for k, v in sorted(c.items(), reverse=True)]
        return f"{len(self.ops)} gates [{', '.join(parts) or 'empty'}], depth {self.depth()}"


def make_circuit(num_qubits: int, ops: Iterable[Gate]) -> Circuit:
    return Circuit(num_qubits, list(ops))
