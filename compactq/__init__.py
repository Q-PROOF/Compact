"""Q-PROOF Compact — the verified quantum circuit optimizer.

Smaller circuits, proven.

Smaller, shallower quantum circuits, exact by construction and verified
before return.  Works standalone or alongside Qiskit.

    >>> import compactq
    >>> from compactq.benchmarks import qft
    >>> opt = compactq.optimize(qft(4))
    >>> opt.two_qubit_count() <= qft(4).two_qubit_count()
    True

The public names load lazily (PEP 562) so `import compactq` costs one
small module instead of the whole stack — CLI and agent callers get
interactive latency; explicit imports (`from compactq import optimize`,
`import compactq.noise`) work exactly as before.
"""
from .circuit import Circuit, Gate

__version__ = "0.2.3"

_EXPORTS = {
    "optimize": ".optimize",
    "optimize_deep": ".optimize",
    "optimize_search": ".search",
    "is_clifford": ".stabilizer",
    "clifford_equal": ".stabilizer",
    "approximate": ".approximate",
    "Target": ".target",
    "optimize_for": ".target",
    "approximate_for_target": ".target",
    "from_qasm3": ".io_qasm",
    "from_qasm": ".io_qasm",
    "to_qasm": ".io_qasm",
    "to_qasm3": ".io_qasm3",
    "expand_mcx": ".mcx",
    "expand_mcp": ".mcx",
    "template_pass": ".templates",
    "param": ".symbolic",
    "bind": ".symbolic",
    "structure_optimize": ".symbolic",
    "optimize_large": ".verify_large",
    "states_agree": ".verify_large",
    "verify": ".verify",
    "VERIFICATION_TIERS": ".verify",
    "build_certificate": ".certificate",
    "verify_compositional": ".compositional",
    "verify_segmented": ".compositional",
    "coupling_preset": ".topology",
    "UnsupportedCircuitError": ".errors",
    "from_qiskit": ".qiskit_bridge",
    "to_qiskit": ".qiskit_bridge",
    "compactq_pass": ".qiskit_bridge",
    "exact_placement": ".hardware",
    "NoiseModel": ".noise",
    "default_model": ".noise",
    "pauli_twirl": ".suppress",
    "insert_dd": ".suppress",
    "suppress_plan": ".suppress",
    "suppress_execute": ".suppress",
    "expand_for_suppression": ".suppress",
    "mitigate_counts": ".mitigate",
    "mitigate_mle": ".mitigate",
    "simulate_counts": ".simulate",
    "SuppressionReport": ".report",
    "RunTarget": ".adapters",
    "qiskit_runtime": ".adapters",
    "braket_device": ".adapters",
    "stab_sample": ".stabsim",
    "fold_global": ".zne",
    "zne_expectation": ".zne",
    "zne_execute": ".zne",
    "cdr_execute": ".cdr",
    "near_clifford_variants": ".cdr",
    "statevector": ".simulate",
    "exact_probabilities": ".simulate",
    "layer_fidelity": ".metrics",
    "eplg": ".metrics",
    "suppression_metrics": ".metrics",
    "maxcut_qaoa": ".solvers",
    "brute_force_maxcut": ".solvers",
    "shadow_snapshots": ".shadows",
    "shadow_estimate_parity": ".shadows",
    "resource_estimate": ".resources",
    "t_depth": ".resources",
    "rebase_cliffordt": ".resources",
}

__all__ = ["Circuit", "Gate", "benchmarks", "__version__"] + sorted(_EXPORTS)


def __getattr__(name):
    if name in _EXPORTS:
        import importlib
        val = getattr(importlib.import_module(_EXPORTS[name], __name__),
                      name)
        # cache on the package: beats the submodule-attribute binding the
        # import system performs for `import compactq.<module>` and keeps
        # `from compactq import optimize` resolving to the FUNCTION
        globals()[name] = val
        return val
    if name == "benchmarks":
        import importlib
        mod = importlib.import_module(".benchmarks", __name__)
        globals()["benchmarks"] = mod
        return mod
    if name.startswith("to_device"):  # kept for __all__ compatibility
        raise AttributeError(name)
    raise AttributeError(
        f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(list(globals()) + __all__))
