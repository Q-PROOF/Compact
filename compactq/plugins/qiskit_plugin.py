"""Qiskit transpiler plugin pass for compactq.

Wraps the compactq optimizer as a Qiskit TransformationPass:

    from compactq.plugins.qiskit_plugin import CompactqPass
    from qiskit.transpiler import PassManager
    out = PassManager([CompactqPass()]).run(qc)

Requires qiskit (imported lazily); compactq itself stays zero-dependency.
"""
from __future__ import annotations


class CompactqPass:
    """Lazy facade: the real pass class is built by `make_compactq_pass`."""

    _verify = True
    _approximate = False
    _min_fidelity = 0.99


def make_compactq_pass(verify: bool = True, approximate: bool = False,
                   min_fidelity: float = 0.99):
    """Build a Qiskit TransformationPass wrapping the compactq optimizer.

    verify=True enforces the exact whole-circuit proof (<= 6/8 qubits).
    approximate=True switches to the bounded-fidelity approximate mode.
    """
    from qiskit.transpiler import TransformationPass
    from qiskit.converters import circuit_to_dag, dag_to_circuit

    from ..approximate import approximate as _approximate
    from ..search import optimize_search

    class _CompactqPass(TransformationPass):
        def __init__(self):
            super().__init__()
            self._verify = verify
            self._approx = approximate
            self._min_fid = min_fidelity

        def run(self, dag):
            from ..qiskit_bridge import from_qiskit, to_qiskit
            circ = from_qiskit(dag_to_circuit(dag))
            if self._approx:
                out = _approximate(circ, min_fidelity=self._min_fid)
            else:
                out = optimize_search(circ, verify=self._verify)
            return circuit_to_dag(to_qiskit(out))

    return _CompactqPass()


def make_suppression_pass(coupling_map=None, noise_model=None, variants: int = 1,
                          dd_sequence: str = "auto", seed: int = 0):
    """Build a Qiskit TransformationPass applying the suppression
    pipeline circuitally: optimize -> expand -> (layout/route when a
    coupling map is given) -> twirl -> dynamical decoupling.  The pass
    emits variant 0 of `suppress_plan` - exactly equivalent to the
    input (proven before return).

    coupling_map: qiskit CouplingMap or list of edges; noise_model: a
    qiskit BackendV2 (calibrations ingested) or None.
    """
    from qiskit.transpiler import TransformationPass
    from qiskit.converters import circuit_to_dag, dag_to_circuit

    class _SuppressionPass(TransformationPass):
        def __init__(self):
            super().__init__()
            self._coupling = (list(coupling_map.get_edges())
                              if hasattr(coupling_map, "get_edges")
                              else (list(coupling_map) if coupling_map else None))
            self._noise_model = noise_model
            self._variants = variants
            self._dd_sequence = dd_sequence
            self._seed = seed
            self.last_report = None

        def run(self, dag):
            from ..qiskit_bridge import from_qiskit, to_qiskit
            from ..noise import NoiseModel
            from ..suppress import suppress_plan
            circ = from_qiskit(dag_to_circuit(dag))
            if self._noise_model is not None and not isinstance(
                    self._noise_model, NoiseModel):
                noise = NoiseModel.from_qiskit_backend(self._noise_model)
            else:
                noise = self._noise_model
            plan = suppress_plan(circ, noise, seed=self._seed,
                                 variants=self._variants,
                                 coupling=self._coupling,
                                 dd_sequence=self._dd_sequence)
            self.last_report = plan["report"]
            return circuit_to_dag(to_qiskit(plan["variants"][0]))

    return _SuppressionPass()


def suppress_qiskit(qc, coupling_map=None, noise_model=None, variants: int = 1,
                    dd_sequence: str = "auto", seed: int = 0):
    """One-call qiskit-native suppression: returns the suppressed
    QuantumCircuit (exactly equivalent to `qc`)."""
    from qiskit.transpiler import PassManager
    pm = PassManager([make_suppression_pass(
        coupling_map=coupling_map, noise_model=noise_model,
        variants=variants, dd_sequence=dd_sequence, seed=seed)])
    return pm.run(qc)
