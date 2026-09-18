"""Trust artifact: the machine-checkable compilation certificate.

`build_certificate(original, optimized)` returns a JSON-serializable dict
that records WHO compiled, WHAT changed, and the VERIFICATION EVIDENCE —
with sha256 hashes of both circuits' canonical QASM so another machine can
independently re-check the claim:

    compactq verify original.qasm optimized.qasm -o certificate.json
"""
from __future__ import annotations

import hashlib

from .io_qasm import to_qasm

__all__ = ["build_certificate", "qasm_hash"]


def qasm_hash(circ) -> str:
    """sha256 of the circuit's canonical OpenQASM 2.0 serialization."""
    return hashlib.sha256(to_qasm(circ).encode("utf-8")).hexdigest()


def build_certificate(original, optimized, target_name: str | None = None,
                      estimated_fidelity: float | None = None,
                      verification: dict | None = None) -> dict:
    """Build the compilation certificate for an optimized circuit.

    `verification` may be supplied when the caller already ran
    `compactq.verify`; otherwise it is computed here.  `target_name` /
    `estimated_fidelity` optionally record the hardware context.
    """
    import compactq
    v = verification if verification is not None else compactq.verify(
        original, optimized)
    return {
        "compiler": "compactq",
        "version": compactq.__version__,
        "input_hash": qasm_hash(original),
        "output_hash": qasm_hash(optimized),
        "input_qubits": original.num_qubits,
        "optimization": {
            "gates_before": len(original.ops),
            "gates_after": len(optimized.ops),
            "two_qubit_before": original.two_qubit_count(),
            "two_qubit_after": optimized.two_qubit_count(),
            "cx_equivalent_before": original.cx_equivalent_count(),
            "cx_equivalent_after": optimized.cx_equivalent_count(),
            "depth_before": original.depth(),
            "depth_after": optimized.depth(),
        },
        "verification": {
            "equivalent": v.get("equivalent"),
            "tier": v.get("tier"),
            "method": v.get("method"),
            "global_phase_ignored": v.get("global_phase_ignored", True),
        },
        **({"hardware": {"target": target_name,
                         "estimated_fidelity": estimated_fidelity}}
           if (target_name or estimated_fidelity is not None) else {}),
    }
