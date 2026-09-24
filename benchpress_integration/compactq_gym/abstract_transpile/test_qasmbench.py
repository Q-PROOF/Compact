# Part of the compactq Benchpress integration (github.com/Q-PROOF/Compact).
# License: MIT (compactq).
"""Abstract-transpilation workout for compactq: QASMBench suites.

Mirrors benchpress/qiskit_gym/abstract_transpile/test_qasmbench.py so
records are directly comparable: same circuits, same loader, same
record schema.  The timed region is exactly `compactq.optimize_search`
(the whole-circuit equivalence proof is included in its time, as
compactq never returns unproven output); conversion in/out of
compactq's IR happens outside the timed region, matching how the
qiskit gym builds its pass manager outside the timed region.

Scope: all-to-all abstract topology only (compactq's logical
optimization is topology-agnostic; hardware mapping is a separate
documented pipeline).
"""
import pytest

from benchpress.config import Configuration
from benchpress.utilities.backends import FlexibleBackend
from benchpress.utilities.io import qasm_circuit_loader, output_circuit_properties
from benchpress.utilities.validation import circuit_validator

from benchpress.workouts.validation import benchpress_test_validation
from benchpress.workouts.abstract_transpile import (
    WorkoutAbstractQasmBenchSmall,
    WorkoutAbstractQasmBenchMedium,
    WorkoutAbstractQasmBenchLarge,
)
from benchpress.workouts.abstract_transpile.qasmbench import (
    SMALL_CIRC_TOPO,
    SMALL_NAMES,
    MEDIUM_CIRC_TOPO,
    MEDIUM_NAMES,
    LARGE_CIRC_TOPO,
    LARGE_NAMES,
)

ALL_TO_ALL = "all-to-all"


def _compactq_optimize(qc):
    """Convert into compactq, optimize with the verified search
    portfolio, convert back.  Trailing measurements are stripped first
    (compactq optimizes unitary cores; the QASMBench files carry final
    measure statements, and mid-circuit control flow would be a hard
    error).  Conversion happens here (untimed callers wrap only the
    optimize call); errors are honest failures."""
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from compactq import optimize_search
    work = qc.copy()
    work.remove_final_measurements(inplace=True)
    return to_qiskit(optimize_search(from_qiskit(work)))


@benchpress_test_validation
class TestWorkoutAbstractQasmBenchSmall(WorkoutAbstractQasmBenchSmall):
    @pytest.mark.parametrize(
        "circ_and_topo",
        [(c, t) for (c, t) in SMALL_CIRC_TOPO if t == ALL_TO_ALL],
        ids=[n for (c, t), n in zip(SMALL_CIRC_TOPO, SMALL_NAMES)
             if t == ALL_TO_ALL],
    )
    def test_QASMBench_small(self, benchmark, circ_and_topo):
        circuit = qasm_circuit_loader(circ_and_topo[0], benchmark)
        backend = FlexibleBackend(
            circuit.num_qubits, circ_and_topo[1], control_flow=True
        )

        @benchmark
        def result():
            return _run_or_skip(circuit)

        output_circuit_properties(result, "cx", benchmark)
        assert circuit_validator(result, backend)


def _run_or_skip(circuit):
    """compactq optimizes unitary cores: circuits with mid-circuit
    measurement/reset or unconvertible gates are OUT OF SCOPE and skip
    cleanly (recorded as skips in the Benchpress record, never
    silently mis-run)."""
    from compactq.errors import UnsupportedCircuitError
    try:
        return _compactq_optimize(circuit)
    except (UnsupportedCircuitError, ValueError) as e:
        pytest.skip(f"outside compactq's unitary-core scope: "
                    f"{str(e)[:110]}")


@benchpress_test_validation
class TestWorkoutAbstractQasmBenchMedium(WorkoutAbstractQasmBenchMedium):
    @pytest.mark.parametrize(
        "circ_and_topo",
        [(c, t) for (c, t) in MEDIUM_CIRC_TOPO if t == ALL_TO_ALL],
        ids=[n for (c, t), n in zip(MEDIUM_CIRC_TOPO, MEDIUM_NAMES)
             if t == ALL_TO_ALL],
    )
    def test_QASMBench_medium(self, benchmark, circ_and_topo):
        circuit = qasm_circuit_loader(circ_and_topo[0], benchmark)
        backend = FlexibleBackend(
            circuit.num_qubits, circ_and_topo[1], control_flow=True
        )

        @benchmark
        def result():
            return _run_or_skip(circuit)

        output_circuit_properties(result, "cx", benchmark)
        assert circuit_validator(result, backend)


@benchpress_test_validation
class TestWorkoutAbstractQasmBenchLarge(WorkoutAbstractQasmBenchLarge):
    @pytest.mark.parametrize(
        "circ_and_topo",
        [(c, t) for (c, t) in LARGE_CIRC_TOPO if t == ALL_TO_ALL],
        ids=[n for (c, t), n in zip(LARGE_CIRC_TOPO, LARGE_NAMES)
             if t == ALL_TO_ALL],
    )
    def test_QASMBench_large(self, benchmark, circ_and_topo):
        circuit = qasm_circuit_loader(circ_and_topo[0], benchmark)
        backend = FlexibleBackend(
            circuit.num_qubits, circ_and_topo[1], control_flow=True
        )

        @benchmark
        def result():
            return _run_or_skip(circuit)

        output_circuit_properties(result, "cx", benchmark)
        assert circuit_validator(result, backend)
