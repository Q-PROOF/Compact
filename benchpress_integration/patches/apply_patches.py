# Part of the compactq Benchpress integration (github.com/Q-PROOF/Compact).
# License: MIT (compactq).
"""Inject the compactq dispatch branches into a Benchpress checkout.

Benchpress dispatches per-gym output recording and validation on
`Configuration.gym_name`; adding a gym upstream means adding one `elif`
branch to each of two files.  This script applies exactly those two
insertions, idempotently, to a local Benchpress source tree — the same
diff we would contribute upstream.

Usage:  python apply_patches.py <benchpress-checkout-root>
Exits 0 when the tree is patched (or already patched), 1 on error.
"""
from __future__ import annotations

import sys
from pathlib import Path

CIRCUIT_OUTPUT_ELIF = """    elif gym_name == "compactq":
        from benchpress.compactq_gym.utils.io import \\
            compactq_output_circuit_properties

        compactq_output_circuit_properties(circuit, two_qubit_gate, benchmark)
"""

VALIDATION_ELIF = """    elif gym_name in ["compactq"]:
        from benchpress.compactq_gym.utils.validation import \\
            compactq_circuit_validation

        compactq_circuit_validation(circuit, backend)
"""

QASM_LOADER_ELIF = """    if gym_name in ["compactq"]:
        # compactq converts from qiskit circuits (from_qiskit), so the
        # qiskit loader is the right input path for this gym
        from benchpress.qiskit_gym.utils.io import qiskit_qasm_loader

        circuit = qiskit_qasm_loader(qasm_file, benchmark)
    elif gym_name == "qiskit":
"""

INPUT_PROPS_ELIF = """    if gym_name in ["compactq"]:
        from benchpress.qiskit_gym.utils.io import qiskit_input_circuit_properties

        qiskit_input_circuit_properties(circuit, benchmark)
    elif gym_name == "qiskit":
"""


def _patch(path: Path, marker: str, insertion: str) -> bool:
    text = path.read_text(encoding="utf-8")
    if "compactq" in text:
        print(f"  already patched: {path.name}")
        return True
    idx = text.find(marker)
    if idx < 0:
        print(f"  MARKER NOT FOUND in {path}")
        return False
    patched = text[:idx] + insertion + text[idx:]
    path.write_text(patched, encoding="utf-8")
    print(f"  patched: {path.name}")
    return True


def _patch_after_gym_line(path: Path, body: str) -> bool:
    """Insert an early-return compactq branch right after the
    `gym_name = Configuration.gym_name` line of a dispatch function."""
    text = path.read_text(encoding="utf-8")
    if "compactq" in text:
        print(f"  already patched: {path.name}")
        return True
    marker = "    gym_name = Configuration.gym_name\n"
    idx = text.find(marker)
    if idx < 0:
        print(f"  MARKER NOT FOUND in {path.name}")
        return False
    insert_at = idx + len(marker)
    patched = text[:insert_at] + body + text[insert_at:]
    path.write_text(patched, encoding="utf-8")
    print(f"  patched: {path.name}")
    return True


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if isinstance(argv, str):
        argv = [argv]
    args = list(argv)
    if len(args) != 1:
        print(__doc__)
        return 1
    # args[0] is the Benchpress REPO root; the python package lives at
    # <repo>/benchpress/
    root = Path(args[0]) / "benchpress"
    ok = True
    ok &= _patch_after_gym_line(
        root / "utilities" / "io" / "circuit_output.py",
        "    if gym_name == \"compactq\":\n"
        "        from benchpress.compactq_gym.utils.io import compactq_output_circuit_properties\n"
        "\n"
        "        compactq_output_circuit_properties(circuit, two_qubit_gate, benchmark)\n"
        "        return\n")
    ok &= _patch_after_gym_line(
        root / "utilities" / "validation" / "validation.py",
        "    if gym_name == \"compactq\":\n"
        "        from benchpress.compactq_gym.utils.validation import compactq_circuit_validation\n"
        "\n"
        "        compactq_circuit_validation(circuit, backend)\n"
        "        return True\n")
    ok &= _patch_after_gym_line(
        root / "utilities" / "io" / "qasm_loader.py",
        "    if gym_name == \"compactq\":\n"
        "        # compactq converts from qiskit circuits, so the qiskit\n"
        "        # loader is the right input path for this gym\n"
        "        from benchpress.qiskit_gym.utils.io import qiskit_qasm_loader\n"
        "\n"
        "        circuit = qiskit_qasm_loader(qasm_file, benchmark)\n"
        "        return circuit\n")
    ok &= _patch_after_gym_line(
        root / "utilities" / "io" / "circuit_input.py",
        "    if gym_name == \"compactq\":\n"
        "        from benchpress.qiskit_gym.utils.io import qiskit_input_circuit_properties\n"
        "\n"
        "        qiskit_input_circuit_properties(circuit, benchmark)\n"
        "        return\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
