"""Qiskit transpiler plugin pass for compactq.

Wraps the compactq optimizer as a Qiskit TransformationPass:

    from compactq.plugins.qiskit_plugin import CompactqPass
    from qiskit.transpiler import PassManager
    out = PassManager([CompactqPass()]).run(qc)

The proof is preserved as inspectable pipeline output, never silently
dropped: after the pass runs, ``PassManager.property_set["compactq_
proof"]`` and the output circuit's ``metadata["compactq"]`` carry the
verification verdict (tier, method, wall time).

Requires qiskit (imported lazily); compactq itself stays zero-dependency.
"""
from __future__ import annotations


class CompactqPass:
    """compactq as a Qiskit pass.

    Instantiating ``CompactqPass(...)`` returns a real
    ``qiskit.TransformationPass`` (qiskit is imported lazily so the core
    package stays zero-dependency).  Keyword arguments:

      optimization_level  None -> full search portfolio (the same
                          candidates as ``optimize_search``); <= 0 -> the
                          plain ``optimize`` pipeline.
      verify=True         enforce the exact whole-circuit proof before
                          returning (default; False is NOT approximate —
                          it skips only the internal re-check, every
                          rewrite is exact either way).
      approximate=False   switch to the bounded-fidelity approximate mode
                          (requires min_fidelity).
      proof_metadata=True record the independent verification verdict in
                          ``property_set["compactq_proof"]`` and the
                          output circuit's ``metadata["compactq"]``.
    """

    _verify = True
    _approximate = False
    _min_fidelity = 0.99

    def __new__(cls, optimization_level: int | None = None,
                verify: bool = True, approximate: bool = False,
                min_fidelity: float = 0.99, proof_metadata: bool = True):
        return make_compactq_pass(
            verify=verify, approximate=approximate,
            min_fidelity=min_fidelity,
            optimization_level=optimization_level,
            proof_metadata=proof_metadata)


def make_compactq_pass(verify: bool = True, approximate: bool = False,
                   min_fidelity: float = 0.99,
                   optimization_level: int | None = None,
                   proof_metadata: bool = True):
    """Build a Qiskit TransformationPass wrapping the compactq optimizer.

    verify=True enforces the exact whole-circuit proof (<= 6/8 qubits).
    approximate=True switches to the bounded-fidelity approximate mode.
    With proof_metadata=True the pass writes the verification verdict
    (``compactq.verify`` output) to ``self.property_set["compactq_proof"]``
    and to the output circuit's ``metadata["compactq"]``; the pass
    instance also keeps it as ``self.last_proof``.
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
            self._opt_level = optimization_level
            self._proof_metadata = proof_metadata
            self.last_proof = None

        def run(self, dag):
            import compactq
            from ..qiskit_bridge import from_qiskit, to_qiskit
            qc_in = dag_to_circuit(dag)
            circ = from_qiskit(qc_in)
            if self._approx:
                out = _approximate(circ, min_fidelity=self._min_fid)
            elif self._opt_level is not None and self._opt_level <= 0:
                from ..optimize import optimize
                out = optimize(circ)
            else:
                out = optimize_search(circ, verify=self._verify)
            out_qc = to_qiskit(out)
            if self._proof_metadata:
                proof = dict(compactq.verify(circ, out))
                proof["approximate"] = bool(self._approx)
                self.last_proof = proof
                ps = getattr(self, "property_set", None)
                if ps is not None:
                    ps["compactq_proof"] = proof
                meta = dict(out_qc.metadata or {})
                meta["compactq"] = proof
                out_qc.metadata = meta
            return circuit_to_dag(out_qc)

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
