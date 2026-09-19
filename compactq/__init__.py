"""Q-PROOF Compact — the verified quantum circuit optimizer.

Smaller circuits, proven.

Smaller, shallower quantum circuits, exact by construction and verified
before return.  Works standalone or alongside Qiskit.

    >>> import compactq
    >>> from compactq.benchmarks import qft
    >>> opt = compactq.optimize(qft(4))
    >>> opt.two_qubit_count() <= qft(4).two_qubit_count()
    True
"""
from .circuit import Circuit, Gate
from .optimize import optimize, optimize_deep
from .search import optimize_search
from .approximate import approximate
from .target import Target, optimize_for, approximate_for_target
from .mcx import expand_mcx, expand_mcp
from .templates import template_pass
from .stabilizer import clifford_equal, is_clifford
from .io_qasm import from_qasm, from_qasm3, to_qasm
from .errors import UnsupportedCircuitError
from .symbolic import param, bind, structure_optimize
from .verify_large import optimize_large, states_agree
from .verify import verify, VERIFICATION_TIERS
from .certificate import build_certificate
from .compositional import verify_compositional, verify_segmented
from .topology import coupling_preset
from .io_qasm3 import to_qasm3
from .qiskit_bridge import from_qiskit, to_qiskit, compactq_pass
from .hardware import exact_placement
from .noise import NoiseModel, default_model
from .suppress import (pauli_twirl, insert_dd, suppress_plan,
                       suppress_execute, expand_for_suppression)
from .mitigate import mitigate_counts, mitigate_mle
from .adapters import RunTarget, qiskit_runtime, braket_device
from .report import SuppressionReport
from .simulate import simulate_counts
from .stabsim import stab_sample
from .zne import fold_global, zne_expectation, zne_execute
from .cdr import cdr_execute, near_clifford_variants
from .metrics import layer_fidelity, eplg, suppression_metrics
from .shadows import shadow_snapshots, shadow_estimate_parity
from .resources import resource_estimate, t_depth, rebase_cliffordt
from .solvers import maxcut_qaoa, brute_force_maxcut
from .simulate import statevector, exact_probabilities
from . import benchmarks

__version__ = "0.2.0"
__all__ = ["Circuit", "Gate", "optimize", "optimize_deep", "optimize_search", "is_clifford", "clifford_equal", "approximate", "Target", "optimize_for", "approximate_for_target", "from_qasm3", "expand_mcx", "expand_mcp", "template_pass", "from_qasm", "to_qasm",
           "to_qasm3", "benchmarks", "__version__", "param", "bind",
            "structure_optimize", "optimize_large", "states_agree", "verify",
            "VERIFICATION_TIERS", "build_certificate", "verify_compositional",
            "verify_segmented", "coupling_preset",
            "UnsupportedCircuitError",
            "NoiseModel", "default_model", "pauli_twirl", "insert_dd",
            "suppress_plan", "suppress_execute", "expand_for_suppression",
            "mitigate_counts", "mitigate_mle", "simulate_counts",
            "SuppressionReport", "to_device",
            "RunTarget", "qiskit_runtime", "braket_device", "stab_sample",
            "fold_global", "zne_expectation", "zne_execute",
            "cdr_execute", "near_clifford_variants", "statevector",
            "exact_probabilities", "layer_fidelity", "eplg",
            "suppression_metrics", "maxcut_qaoa", "brute_force_maxcut",
            "shadow_snapshots", "shadow_estimate_parity",
            "resource_estimate", "t_depth", "rebase_cliffordt",
            "exact_placement"]
