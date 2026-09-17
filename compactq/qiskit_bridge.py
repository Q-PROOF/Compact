"""Optional Qiskit interoperability (requires `pip install qiskit`).

Provides loss-free circuit conversion and a transpiler-compatible pass so
compactq can run inside an existing Qiskit workflow:

    from compactq.qiskit_bridge import compactq_pass, to_qiskit, from_qiskit
    qc2 = compactq_pass(qc)                      # optimize a QuantumCircuit
"""
from __future__ import annotations

from .circuit import Circuit, Gate, ONE_Q_GATES, TWO_Q_GATES, MULTI_Q_GATES
from .optimize import optimize
from .search import optimize_search

_CONVERT_1Q = {"h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p", "sx"}

# Names the compactq IR can represent natively; anything else a QuantumCircuit
# carries (mcphase/ccx/ecr/rxx/iswap/cu3/...) is decomposed at this boundary
# into the standard-gate basis below via qiskit's own equivalence library, so
# conversion stays loss-free instead of crashing on gate construction.
_COMPACTQ_BASIS = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz",
               "p", "sx", "sxdg", "u", "cx", "cz", "swap", "cp", "ecr", "iswap"]
_SKIP_NAMES = {"barrier"}
_NON_UNITARY_NAMES = {"measure", "reset", "if_else", "while_loop",
                      "for_loop", "switch_case", "store", "set_register",
                      "expression", "linear_function"}


def to_qiskit(circ: Circuit):
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(circ.num_qubits)
    for g in circ.ops:
        params = list(g.params)
        if len(g.qubits) == 2 and g.name == "swap":
            qc.swap(g.qubits[0], g.qubits[1])
        elif g.name == "u3":  # QuantumCircuit.u3 was removed in Qiskit 2.x
            qc.u(*params, *g.qubits)
        else:
            getattr(qc, g.name)(*params, *g.qubits)
    return qc


def from_qiskit(qc) -> Circuit:
    """Convert a QuantumCircuit into compactq's IR, loss-free.

    Raises UnsupportedCircuitError when the circuit contains measurement,
    reset, or classical control flow: those operations change program
    semantics, and silently dropping them would return a *different*
    program.  Strip them first, e.g.
    ``qc.remove_final_measurements(inplace=True)`` for trailing
    measurements.  Barriers pass through as the no-ops they are.
    """
    from .errors import UnsupportedCircuitError
    non_unitary = sorted({inst.operation.name for inst in qc.data}
                         & _NON_UNITARY_NAMES)
    if non_unitary:
        raise UnsupportedCircuitError(
            "circuit contains non-unitary operation(s): "
            + ", ".join(non_unitary)
            + ". compactq optimizes unitary circuits; strip them first "
            "(e.g. qc.remove_final_measurements(inplace=True) for "
            "trailing measurements).")
    if not (set(names_safe(qc)) <= ONE_Q_GATES | TWO_Q_GATES | MULTI_Q_GATES):
        from qiskit import transpile
        try:
            qc = transpile(qc, basis_gates=_COMPACTQ_BASIS,
                           optimization_level=0, seed_transpiler=42)
        except Exception:
            # exotic gates that qiskit itself cannot rebase: fall through
            # to gate construction, which will raise a precise error
            pass
    circ = Circuit(qc.num_qubits, [])
    for inst in qc.data:
        name = inst.operation.name
        if name == "barrier":
            continue
        qubits = tuple(qc.find_bit(q).index for q in inst.qubits)
        params = []
        for p in inst.operation.params:
            try:
                params.append(float(p))
            except (TypeError, ValueError):
                pass
        circ.append(Gate(name, tuple(params), qubits))
    return circ


def names_safe(qc):
    return {inst.operation.name for inst in qc.data} - _SKIP_NAMES


def compactq_pass(qc, optimization_level: int | None = None):
    """Optimize a Qiskit QuantumCircuit with compactq and return a QuantumCircuit.

    Uses the multi-pipeline search (verified for circuits within the
    prover's qubit limit); ``optimization_level`` is accepted for
    transpiler-plugin compatibility and selects plain optimize at level 0.
    """
    circ = from_qiskit(qc)
    if optimization_level is not None and optimization_level <= 0:
        opt = optimize(circ)
    else:
        opt = optimize_search(circ)
    return to_qiskit(opt)
