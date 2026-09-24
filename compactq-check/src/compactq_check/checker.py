"""Certificate validation: re-derives every claim from first principles.

Validates a Compact compilation certificate against the two provided
QASM circuits:
  1. schema check (version, composition, witness kind),
  2. hash pinning (sha256 of the provided QASM must match the certificate),
  3. independent re-computation:
       whole_circuit_v1      per witness kind:
         unitary           — dense rebuild of both unitaries (≤ 8 qubits)
         stabilizer        — fresh Clifford canonical forms, any width
         phase_polynomial  — fresh parity→angle tables, any width
         dd                — fresh decision-diagram overlap (≤ 32 qubits,
                             node-budgeted; the checker declines, never
                             guesses, when its budget is exhausted)
       disjoint_blocks_v1    — the certificate's block QASMs are re-parsed
                             and re-proven dense pair-by-pair; the block
                             structure (disjoint, covering, identical
                             between the two circuits) is re-checked here
       sequential_segments_v1 — the certificate's segments must tile both
                             provided circuits op-for-op, then every
                             segment pair is re-proven dense on its
                             active qubits

Verdicts: VALID / INVALID (well-formed but false or unprovable) —
MALFORMED is raised as ValueError by the loader for unparseable inputs.
"""
from __future__ import annotations

import hashlib

from .ops import (active_components, build_unitary, build_unitary_active,
                  clifford_equal, dd_equal, fidelity, is_clifford_ops,
                  phasepoly_canonical, phasepoly_equal)

WITNESS_LIMITS = {"unitary": 8, "stabilizer": None, "phase_polynomial": None,
                  "dd": 32, "composition": None}
COMPOSITIONS = ("whole_circuit_v1", "disjoint_blocks_v1",
                "sequential_segments_v1")
PIECE_LIMIT = 8            # dense re-derivation limit per block/segment


def check_certificate(cert: dict, input_qasm: str, output_qasm: str,
                      tol: float = 1e-7) -> dict:
    problems = []
    if not isinstance(cert, dict):
        return _verdict("MALFORMED", "certificate is not an object")
    if cert.get("cert_version") != 1:
        problems.append(f"cert_version must be 1, got "
                        f"{cert.get('cert_version')!r}")
    composition = cert.get("composition")
    if composition not in COMPOSITIONS:
        problems.append(f"unsupported composition {composition!r}")
    w = cert.get("witness") or {}
    kind = w.get("kind")
    if kind not in WITNESS_LIMITS:
        problems.append(f"unknown witness kind {kind!r}")
    if problems:
        return _verdict("MALFORMED", "; ".join(problems))

    import compactq_check.qasm as q
    n_in, ops_in = q.parse_qasm(input_qasm)
    n_out, ops_out = q.parse_qasm(output_qasm)

    if composition == "disjoint_blocks_v1" and kind == "composition":
        return _check_blocks(cert, input_qasm, output_qasm, q, tol)
    if composition == "sequential_segments_v1" and kind == "composition":
        return _check_segments(cert, input_qasm, output_qasm, q, ops_in,
                               ops_out, tol)
    if kind == "composition":
        return _verdict("MALFORMED",
                        "composition witness requires a compositional "
                        "composition kind")

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


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _pin_full_circuits(cert, input_qasm, output_qasm, n_in):
    """Hash/n_qubit pinning shared by the compositional handlers."""
    if cert.get("input", {}).get("sha256") != _sha(input_qasm):
        return _verdict("INVALID", "input hash mismatch")
    if cert.get("output", {}).get("sha256") != _sha(output_qasm):
        return _verdict("INVALID", "output hash mismatch")
    if n_in != cert.get("n_qubits"):
        return _verdict("INVALID",
                        f"n_qubits {cert.get('n_qubits')} does not match "
                        f"circuits ({n_in})")
    return None


def _check_blocks(cert, input_qasm, output_qasm, q, tol):
    """disjoint_blocks_v1: the checker RE-DERIVES the component
    structure from the provided full circuits itself (union-find over
    gate wires — the certificate's declared pieces are never trusted
    for the proof), requires the same disjoint structure on both sides,
    and re-proves each component pair dense on its own wires."""
    import compactq_check.qasm as q
    n_in, ops_in = q.parse_qasm(input_qasm)
    _n_out, ops_out = q.parse_qasm(output_qasm)
    blocked = _pin_full_circuits(cert, input_qasm, output_qasm, n_in)
    if blocked:
        return blocked
    comps_in = active_components(ops_in)
    comps_out = active_components(ops_out)
    if len(comps_in) < 2 or len(comps_in) != len(comps_out):
        return _verdict("INVALID",
                        f"component structure does not support disjoint-"
                        f"block composition ({len(comps_in)} vs "
                        f"{len(comps_out)} components)")
    for i, ((ws_a, ops_a), (ws_b, ops_b)) in enumerate(
            zip(comps_in, comps_out)):
        if ws_a != ws_b:
            return _verdict("INVALID",
                            f"component {i} wire sets differ between the "
                            f"circuits: {ws_a} vs {ws_b}")
        if len(ws_a) > PIECE_LIMIT:
            return _verdict("INVALID",
                            f"component {i} spans {len(ws_a)} wires — "
                            f"exceeds the {PIECE_LIMIT}-wire dense "
                            f"re-derivation limit")
        ua = build_unitary(ops_a, len(ws_a))
        ub = build_unitary(ops_b, len(ws_b))
        f = fidelity(ua, ub)
        if f < 1 - tol:
            return _verdict("INVALID",
                            f"component {i} not equivalent (fidelity "
                            f"{f:.12f})")
    return {"verdict": "VALID", "equivalent": True,
            "tier_kind": "composition",
            "method": (f"disjoint_blocks ({len(comps_in)} components, "
                       f"dense per component, re-derived independently)"),
            "global_phase_ignored": True}


def _check_segments(cert, input_qasm, output_qasm, q, ops_in, ops_out, tol):
    """sequential_segments_v1: the certificate's segments must tile both
    provided circuits op-for-op, then every segment pair is re-proven
    dense on its active wires."""
    n_in, _ = q.parse_qasm(input_qasm)
    blocked = _pin_full_circuits(cert, input_qasm, output_qasm, n_in)
    if blocked:
        return blocked
    segs = cert.get("segments") or {}
    orig, out = segs.get("original"), segs.get("output")
    if not (isinstance(orig, list) and isinstance(out, list)):
        return _verdict("INVALID", "segments.original/segments.output missing")
    if len(orig) < 2 or len(orig) != len(out):
        return _verdict("INVALID",
                        "segments must list >= 2 pieces on both sides, "
                        "with equal counts")

    def parse_tiling(side, full_ops):
        pieces = []
        tiled = []
        for rec in side:
            if not isinstance(rec, dict) or "qasm" not in rec:
                return None, _verdict("INVALID", "malformed segment record")
            if rec.get("sha256") != _sha(rec["qasm"]):
                return None, _verdict("INVALID", "segment hash mismatch")
            _n, ops_p = q.parse_qasm(rec["qasm"])
            if not ops_p:
                return _verdict("INVALID", "empty segment")
            pieces.append(ops_p)
            tiled.extend(ops_p)
        if tiled != full_ops:
            return None, _verdict(
                "INVALID",
                "segments do not tile the provided circuit op-for-op")
        return pieces, None

    pieces_in, err = parse_tiling(orig, ops_in)
    if err:
        return err
    pieces_out, err = parse_tiling(out, ops_out)
    if err:
        return err
    for i, (ops_a, ops_b) in enumerate(zip(pieces_in, pieces_out)):
        wa, ua = build_unitary_active(ops_a, n_in)
        wb, ub = build_unitary_active(ops_b, n_in)
        if wa != wb:
            return _verdict("INVALID",
                            f"segment {i} active wires differ between the "
                            f"circuits: {wa} vs {wb}")
        if len(wa) > PIECE_LIMIT:
            return _verdict("INVALID",
                            f"segment {i} spans {len(wa)} active wires — "
                            f"exceeds the {PIECE_LIMIT}-wire dense "
                            f"re-derivation limit")
        f = fidelity(ua, ub)
        if f < 1 - tol:
            return _verdict("INVALID",
                            f"segment {i} not equivalent (fidelity "
                            f"{f:.12f})")
    return {"verdict": "VALID", "equivalent": True,
            "tier_kind": "composition",
            "method": (f"sequential_segments ({len(orig)} pieces, dense "
                       f"on active wires)"),
            "global_phase_ignored": True}
