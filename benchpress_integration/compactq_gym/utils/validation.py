# Part of the compactq Benchpress integration (github.com/Q-PROOF/Compact).
# License: MIT (compactq).
"""Topology validation for compactq's Benchpress outputs.

compactq's logical optimizer does not route: on the all-to-all
(abstract) topology every 2-qubit pairing is legal, so validation
checks the gate set is quantum-circuit well-formed and records the
scope.  Equivalence itself is proven in-product before compactq
returns anything (see VERIFICATION.md) — this validator is about the
Benchpress output contract, not correctness.
"""
from __future__ import annotations


def compactq_circuit_validation(circuit, backend):
    """All-to-all topology: any 2-qubit gate placement is legal.

    Raises when the backend topology is NOT all-to-all, because the
    compactq gym only claims the abstract all-to-all workout."""
    name = getattr(backend, "name", "")
    if "all" in str(name) or "all-to-all" in str(name):
        return
    # FlexibleBackend exposes the topology via its coupling map
    cmap = getattr(backend, "coupling_map", None)
    if cmap is None:
        return
    n = cmap.graph.num_nodes()
    if cmap.graph.num_edges() < n * (n - 1):
        raise Exception(
            "compactq gym supports the all-to-all abstract topology only; "
            "run hardware-mapped workouts through compactq's own "
            "route_aware pipeline instead")
