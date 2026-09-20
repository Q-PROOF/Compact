"""Certificate validation: re-derives every claim from first principles.

Validates a Compact compilation certificate against the two provided
QASM circuits:
  1. schema check (version, composition, witness kind),
  2. hash pinning (sha256 of the provided QASM must match the certificate),
  3. independent re-computation per witness kind:
       unitary           — dense rebuild of both unitaries (≤ 8 qubits)
       stabilizer        — fresh Clifford canonical forms, any width
       phase_polynomial  — fresh parity→angle tables, any width
       dd                — fresh decision-diagram overlap (≤ 32 qubits,
                           node-budgeted; the checker declines, never
                           guesses, when its budget is exhausted)

Verdicts: VALID / INVALID (well-formed but false or unprovable) —
MALFORMED is raised as ValueError by the loader for unparseable inputs.
"""
from __future__ import annotations

import hashlib

from .ops import (build_unitary, clifford_equal, dd_equal, fidelity,
                  is_clifford_ops, phasepoly_canonical, phasepoly_equal)

WITNESS_LIMITS = {"unitary": 8, "stabilizer": None, "phase_polynomial": None,
                  "dd": 32}


def check_certificate(cert: dict, input_qasm: str, output_qasm: str,
                      tol: float = 1e-7) -> dict:
    problems = []
    if not isinstance(cert, dict):
        return _verdict("MALFORMED", "certificate is not an object")
    if cert.get("cert_version") != 1:
        problems.append(f"cert_version must be 1, got "
                        f"{cert.get('cert_version')!r}")
    if cert.get("composition") != "whole_circuit_v1":
        problems.append(f"unsupported composition "
                        f"{cert.get('composition')!r}")
    w = cert.get("witness") or {}
    kind = w.get("kind")
    if kind not in WITNESS_LIMITS:
        problems.append(f"unknown witness kind {kind!r}")
    if problems:
        return _verdict("MALFORMED", "; ".join(problems))

    import compactq_check.qasm as q
    n_in, ops_in = q.parse_qasm(input_qasm)
    n_out, ops_out = q.parse_qasm(output_qasm)

    if n_in != n_out:
        return _verdict("INVALID", f"qubit count mismatch: {n_in} vs {n_out}")
    if n_in != cert.get("n_qubits"):
        problems.append(f"n_qubits {cert.get('n_qubits')} does not match "
                        f"circuits ({n_in})")

    h_in = hashlib.sha256(input_qasm.encode("utf-8")).hexdigest()
    h_out = hashlib.sha256(output_qasm.encode("utf-8")).hexdigest()
    if cert.get("input", {}).get("sha256") != h_in:
        problems.append("input hash mismatch")
    if cert.get("output", {}).get("sha256") != h_out:
        problems.append("output hash mismatch")
    if problems:
        return _verdict("INVALID", "; ".join(problems))

    limit = WITNESS_LIMITS[kind]
    if limit is not None and n_in > limit:
        return _verdict("INVALID",
                        f"{kind} witness violates the {limit}-qubit limit")

    equivalent = None
    method = kind
    if kind == "unitary":
        ua = build_unitary(ops_in, n_in)
        ub = build_unitary(ops_out, n_in)
        f = fidelity(ua, ub)
        equivalent = f >= 1 - tol
        method = f"full_unitary (fidelity {f:.12f})"
    elif kind == "stabilizer":
        equivalent = clifford_equal(ops_in, ops_out, n_in)
        method = "clifford_tableau_equality"
    elif kind == "phase_polynomial":
        if not (is_phasepoly_ops(ops_in) and is_phasepoly_ops(ops_out)):
            return _verdict("INVALID",
                            "phase_polynomial witness on non-phase-"
                            "polynomial circuits")
        equivalent = phasepoly_equal(ops_in, ops_out, n_in, tol)
        method = "phase_polynomial_table_equality"
    elif kind == "dd":
        eq = dd_equal(ops_in, ops_out, n_in, tol)
        if eq is None:
            return _verdict("INVALID",
                            "dd witness exceeds the checker's own node "
                            "budget — cannot be independently re-derived")
        equivalent = eq
        method = "decision_diagram_overlap"

    if equivalent is True:
        return {"verdict": "VALID", "equivalent": True, "tier_kind": kind,
                "method": method,
                "global_phase_ignored": True}
    return {"verdict": "INVALID", "equivalent": bool(equivalent),
            "tier_kind": kind, "method": method,
            "global_phase_ignored": True}


def _verdict(state: str, detail: str) -> dict:
    return {"verdict": state, "detail": detail,
            "equivalent": None, "tier_kind": None, "method": None,
            "global_phase_ignored": True}
