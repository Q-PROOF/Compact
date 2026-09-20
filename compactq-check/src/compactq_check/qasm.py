"""Minimal, strict OpenQASM 2.0 reader for the checker's gate subset.

Deliberately independent of compactq's parser.  Supports: OPENQASM 2.0
header, qreg, barrier (ignored), and the builtin one/two/three-qubit
gates used by Compact's emitter.  Custom gate definitions are NOT
supported and are reported as unsupported (the certificate spec restricts
certificates to the builtin subset).
"""
from __future__ import annotations

import re

GATE_RE = re.compile(
    r"^\s*([a-zA-Z][a-zA-Z0-9_]*)\s*(?:\(([^)]*)\))?\s*"
    r"((?:q\[\d+\]|c\[\d+\])(?:\s*,\s*(?:q\[\d+\]|c\[\d+\]))*)\s*;\s*$")
QREG_RE = re.compile(r"^\s*qreg\s+([a-zA-Z][a-zA-Z0-9_]*)\s*\[\s*(\d+)\s*\]\s*;\s*$")
CREG_RE = re.compile(r"^\s*creg\s+", re.IGNORECASE)
IDENT_RE = re.compile(r"^\s*([a-zA-Z][a-zA-Z0-9_]*)\s*(?:\(([^)]*)\))?\s*"
                      r"((?:[a-zA-Z][a-zA-Z0-9_]*\[\d+\])(?:\s*,\s*"
                      r"[a-zA-Z][a-zA-Z0-9_]*\[\d+\])*)\s*;\s*$")

SUPPORTED = {
    "id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg",
    "rx", "ry", "rz", "p", "u1", "u2", "u3", "u",
    "cx", "cy", "cz", "swap", "ccx", "cp", "cu1",
}


def parse_qasm(text: str):
    """Parse a gate-subset OpenQASM 2.0 file.

    Returns (n_qubits, [(name, [params], [wires])]).
    Raises ValueError on unsupported/malformed content.
    """
    n = None
    ops = []
    regname = None
    for raw in text.splitlines():
        line = raw.split("//")[0].strip()
        if not line:
            continue
        if line.startswith("OPENQASM"):
            if not line.startswith("OPENQASM 2"):
                raise ValueError("only OPENQASM 2.0 is supported")
            continue
        if line.startswith("include"):
            continue
        if CREG_RE.match(line):
            continue  # classical registers carry no unitary semantics
        m = QREG_RE.match(line)
        if m:
            if n is not None:
                raise ValueError("multiple qreg declarations unsupported")
            regname, n = m.group(1), int(m.group(2))
            continue
        if line.startswith(("creg", "measure", "reset", "if")):
            raise ValueError(f"non-unitary statement unsupported: {line}")
        m = IDENT_RE.match(line)
        if not m:
            raise ValueError(f"unparseable line: {raw.strip()!r}")
        name = m.group(1)
        if name == "barrier":
            continue
        if name not in SUPPORTED:
            raise ValueError(f"unsupported gate: {name}")
        params = []
        if m.group(2):
            for p in m.group(2).split(","):
                p = p.strip()
                try:
                    params.append(float(p))
                except ValueError:
                    raise ValueError(f"non-numeric parameter {p!r}")
        qubits = []
        for q in re.findall(r"[\w]+\[(\d+)\]", m.group(3)):
            idx = int(q)
            if regname is not None and idx >= n:
                raise ValueError("qubit index out of range")
            qubits.append(idx)
        expected = {"cx": 2, "cy": 2, "cz": 2, "swap": 2, "cp": 2,
                    "cu1": 2, "ccx": 3}.get(name, 1)
        if len(qubits) != expected:
            raise ValueError(f"{name} expects {expected} qubits, "
                             f"got {len(qubits)}")
        ops.append((name, params, qubits))
    if n is None:
        raise ValueError("no qreg declaration found")
    return n, ops
