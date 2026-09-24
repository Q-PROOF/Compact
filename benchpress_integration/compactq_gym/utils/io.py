# Part of the compactq Benchpress integration (github.com/Q-PROOF/Compact).
# License: MIT (compactq).
"""compactq-side circuit statistics for Benchpress records.

Key names match the qiskit gym's schema exactly (output_num_qubits /
output_circuit_operations / output_gate_count_2q / output_depth_2q) so
records from the two gyms are directly comparable; the compactq-specific
aliases are written alongside for downstream convenience.
"""
from __future__ import annotations


def compactq_output_circuit_properties(circuit, two_qubit_gate, benchmark):
    """Record the optimized circuit's statistics.

    `circuit` is a qiskit QuantumCircuit (compactq's output converted at
    the boundary, outside the timed region)."""
    ops = {}
    two_qubit_gate_count = 0
    for inst in circuit.data:
        nm = inst.operation.name
        ops[nm] = ops.get(nm, 0) + 1
        if len(inst.qubits) == 2 and nm == two_qubit_gate:
            two_qubit_gate_count += 1
    benchmark.extra_info["output_num_qubits"] = circuit.num_qubits
    benchmark.extra_info["output_circuit_operations"] = ops
    benchmark.extra_info["output_gate_count_2q"] = two_qubit_gate_count
    benchmark.extra_info["output_depth_2q"] = circuit.depth(
        filter_function=lambda x: len(x.qubits) == 2
    )
    # compactq-specific aliases (superset of the qiskit gym schema)
    benchmark.extra_info["output_gate_count"] = circuit.size()
    benchmark.extra_info["output_depth"] = circuit.depth()
    benchmark.extra_info["output_2q_gate_count"] = two_qubit_gate_count
