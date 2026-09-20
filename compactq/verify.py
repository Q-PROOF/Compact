"""Public verification API — the evidence layer of verified compilation.

`verify(original, optimized)` independently checks that two circuits
implement the same unitary (up to global phase) and reports WHICH grade of
evidence backs the answer:

  tier 3  clifford_tableau     both circuits are Clifford: exact algebraic
                               tableau proof, any qubit count
  tier 3  phase_polynomial     both circuits are CX+diagonal: exact GF(2)
                               parity-table proof, any qubit count
  tier 2  full_unitary         dense |Tr(U+V)|/d comparison within the
                               prover's qubit limit (8 with the native
                               kernels, 6 without)
  tier 2  dd_full_unitary      decision-diagram |Tr(U+V)|/d beyond the
                               dense ceiling (structured circuits to 25+
                               qubits; node-budgeted, declines on blowup)
  tier 1  randomized_sampling  K random product states (needs numpy;
                               probabilistically exact to 1e-8)
  tier 0  none                 prover unavailable (e.g. numpy missing on a
                               >limit circuit) — `equivalent` is None

Tier 4 (compositional certificates — disjoint-block decomposition, and
sequential segments via `verify_segmented`) ships today; tier 5 (formal
proof) is roadmap.  The tier field is part of the stable API so callers
can require a minimum grade.

Guarantee vocabulary (used by VERIFICATION_TIERS and the CLI): tier 2/3
are **exact-proven** (numerical full-unitary / algebraic tableau
certificates); tier 1 is **verified-randomized** (probabilistically
exact K-state sampling); tier 0 is **unverified**.  "Exact-proven" here
means a machine-checked certificate of unitary equivalence up to global
phase — it is NOT a formal proof-assistant certificate (that is tier 5).
"""
from __future__ import annotations

import time

from .circuit import Circuit

__all__ = ["verify", "VERIFICATION_TIERS", "proof_status_tier"]

VERIFICATION_TIERS = {
    0: "unverified",
    1: "verified-randomized",
    2: "exact-proven (numerical full-unitary)",
    3: "exact-proven (algebraic tableau)",
    4: "exact-proven (compositional certificates)",
    5: "formal-proof (roadmap)",
}

# CLI proof-status strings -> verification tier (see compactq.cli)
_STATUS_TIER = {
    "exact-unitary": 2,
    "exact-proven (algebraic)": 3,
    "exact-proven (decision-diagram)": 2,
    "compositional": 4,
    "randomized-exact": 1,
    "approximate": 0,
    "unverified": 0,
    "verification-rejected (original returned)": 0,
}


def proof_status_tier(status: str) -> int:
    """Map a CLI proof-status string to its verification tier (0-4)."""
    return _STATUS_TIER.get(status, 0)


def verify(original: Circuit, optimized: Circuit) -> dict:
    """Independently verify that `optimized` implements the same unitary as
    `original`, and return the evidence.

    Returns a dict:
      equivalent          True / False / None (None = could not decide)
      tier                0-4 today (5 = formal proof, roadmap; see
                          VERIFICATION_TIERS).  Tier 4 fires only when the
                          circuit structure satisfies the compositional
                          strategy (disjoint blocks); sequential-segment
                          proofs are available via
                          compactq.compositional.verify_segmented.
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

    # tier 3: CX+diagonal fragment -> exact GF(2) parity-table proof,
    # any width (decisive in both directions when both are in fragment)
    from .phasepoly import phasepoly_equal
    try:
        pp = phasepoly_equal(original, optimized)
        if pp is not None:
            return pack(pp, 3, "phase_polynomial")
    except Exception:
        pass  # fall through on numerical doubt

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

    # tier 2: decision-diagram proof beyond the dense ceiling — structured
    # circuits prove exactly to 25+ qubits; the node budget makes a
    # blowup a decline (fall through), never a wrong answer
    from .dd import DD_MAX_QUBITS, check_equivalent_dd, DDOverflow
    if original.num_qubits <= DD_MAX_QUBITS:
        try:
            return pack(bool(check_equivalent_dd(original, optimized)),
                        2, "dd_full_unitary")
        except Exception:
            pass  # budget exhausted / unsupported gate: fall through

    # tier 1: randomized K-state sampling (numpy optional)
    try:
        from .verify_large import states_agree
        return pack(bool(states_agree(original, optimized)),
                    1, "randomized_sampling")
    except ImportError:
        return pack(None, 0, "prover_unavailable")
