"""External tool interfaces: stim export + optional QCEC referee.

Stim exporter: Clifford-only circuits can be exported to Stim's text
format for high-performance QEC simulation.  Validation requires the
optional `stim` package.

QCEC referee: when mqt.qcec is installed, provides independent
equivalence checking as an additional referee alongside compactq's
built-in proof net.  Optional - never required.
"""
from __future__ import annotations


def to_stim(circ) -> str:
    """Export a Clifford-only circuit to Stim's text format.
    Raises ValueError for non-Clifford gates."""
    lines = []
    n = circ.num_qubits
    for g in circ.ops:
        nm = g.name
        qs = " ".join(str(q) for q in g.qubits)
        if nm == "cx":
            lines.append(f"CX {qs}")
        elif nm == "cz":
            lines.append(f"CZ {qs}")
        elif nm in ("x", "y", "z", "h", "s", "sdg"):
            lines.append(f"{nm.upper()} {qs}")
        elif nm == "swap":
            lines.append(f"SWAP {qs}")
        else:
            raise ValueError(f"gate {nm} is not exportable to Stim")
    lines.append(f"M {n}")
    return "\n".join(lines)


def qcec_referee(circ_a, circ_b) -> bool | None:
    """Independent equivalence check via MQT QCEC (optional import).
    Returns True/False, or None when QCEC is not installed."""
    try:
        from mqt import qcec
    except Exception:
        return None
    from .qiskit_bridge import to_qiskit
    qc_a = to_qiskit(circ_a)
    qc_b = to_qiskit(circ_b)
    result = qcec.verify(qc_a, qc_b)
    if hasattr(result, 'is_equivalent'):
        return result.is_equivalent
    return None
