"""Public verification API — the evidence layer of verified compilation.

`verify(original, optimized)` independently checks that two circuits
implement the same unitary (up to global phase) and reports WHICH grade of
evidence backs the answer:

  tier 3  clifford_tableau     both circuits are Clifford: exact algebraic
                               tableau proof, any qubit count
  tier 2  full_unitary         dense |Tr(U+V)|/d comparison within the
                               prover's qubit limit (8 with the native
                               kernels, 6 without)
  tier 1  randomized_sampling  K random product states (needs numpy;
                               probabilistically exact to 1e-8)
  tier 0  none                 prover unavailable (e.g. numpy missing on a
                               >limit circuit) — `equivalent` is None

Tiers 4 (compositional certificates) and 5 (formal proof) are roadmap; the
tier field is part of the stable API so callers can require a minimum grade.
"""
from __future__ import annotations

import time

from .circuit import Circuit

__all__ = ["verify", "VERIFICATION_TIERS"]

VERIFICATION_TIERS = {
    0: "none",
    1: "randomized_sampling",
    2: "full_unitary",
    3: "clifford_tableau",
    4: "compositional_certificates (roadmap)",
    5: "formal_proof (roadmap)",
}


def verify(original: Circuit, optimized: Circuit) -> dict:
    """Independently verify that `optimized` implements the same unitary as
    `original`, and return the evidence.

    Returns a dict:
      equivalent          True / False / None (None = could not decide)
      tier                0-3 today (4-5 reserved; see VERIFICATION_TIERS)
      method              which prover produced the answer
      global_phase_ignored  always True: equivalence is up to global phase
      runtime_ms          wall time of the check
    """
    t0 = time.perf_counter()

    def pack(equivalent, tier, method, extra=None):
        out = {"equivalent": equivalent, "tier": tier, "method": method,
               "global_phase_ignored": True,
               "runtime_ms": round((time.perf_counter() - t0) * 1000)}
        if extra:
            out.update(extra)
        return out

    if original.num_qubits != optimized.num_qubits:
        return pack(False, 0, "shape_mismatch")

    # tier 3: both Clifford -> exact algebraic tableau proof, any width
    from .stabilizer import is_clifford, clifford_equal
    try:
        if is_clifford(original) and is_clifford(optimized):
            return pack(bool(clifford_equal(original, optimized)),
                        3, "clifford_tableau")
    except Exception:
        pass  # fall through to the numerical provers

    # tier 2: dense unitary comparison within the prover's limit
    from .equivalence import check_equivalent, _MAX_QUBITS
    if original.num_qubits <= _MAX_QUBITS:
        try:
            return pack(bool(check_equivalent(original, optimized)),
                        2, "full_unitary")
        except Exception:
            pass  # numerical doubt: fall through

    # tier 4: compositional proof for block-structured circuits at any
    # width (per-block dense verification over disjoint components — no
    # 2^total unitary is ever built)
    from .compositional import verify_compositional
    comp = verify_compositional(original, optimized)
    if comp is not None:
        extra = {"blocks": comp["blocks"]}
        if not comp["equivalent"]:
            extra["failing_block_qubits"] = comp["failing_block_qubits"]
        return pack(comp["equivalent"], 4, comp["method"], extra)

    # tier 1: randomized K-state sampling (numpy optional)
    try:
        from .verify_large import states_agree
        return pack(bool(states_agree(original, optimized)),
                    1, "randomized_sampling")
    except ImportError:
        return pack(None, 0, "prover_unavailable")
