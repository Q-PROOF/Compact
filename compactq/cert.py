"""Certificates: machine-checkable evidence attached to optimization runs.

Public API:
  optimize_with_certificate(circ) -> CertResult
      .circuit      — the optimized compactq.Circuit
      .certificate  — dict per docs/CERTIFICATE_SPEC.md, or None
      .reason       — why no certificate was issued (when None)
  write_certificate(cert_result, path) — write JSON to a validated path

The certificate's witness is produced by compactq's own provers; the
independent re-check is performed by the separate `compactq-check`
program, which shares no code with this package.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .circuit import Circuit
from .io_qasm import to_qasm

__all__ = ["CertResult", "optimize_with_certificate", "build_certificate",
           "build_certificate_blocks", "build_certificate_segments",
           "write_certificate", "strongest_witness", "CERT_VERSION"]

CERT_VERSION = 1
_COMPOSITION_BLOCK_LIMIT = 8       # per-piece dense prover limit


@dataclass
class CertResult:
    circuit: Circuit
    certificate: Optional[dict]
    reason: str = ""

    def certificate_json(self) -> str:
        if self.certificate is None:
            raise ValueError("no certificate was issued: " + self.reason)
        return json.dumps(self.certificate, indent=2) + "\n"


def write_certificate(result: "CertResult", path: str) -> str:
    """Write the certificate JSON to `path` (expanded + resolved; the
    parent directory must already exist)."""
    out = Path(path).expanduser().resolve()
    if not out.parent.is_dir():
        raise ValueError(f"certificate directory does not exist: "
                         f"{out.parent}")
    out.write_text(result.certificate_json(), encoding="utf-8")
    return str(out)


def _qasm_sha(circ: Circuit) -> str:
    return hashlib.sha256(to_qasm(circ).encode("utf-8")).hexdigest()


def _is_phase_polynomial(circ: Circuit) -> bool:
    """True when every gate is CNOT or diagonal (Z-rotations and Paulis):
    the circuit implements a phase polynomial + linear map."""
    ok = {"cx", "rz", "p", "s", "sdg", "t", "tdg", "z", "id"}
    return all(g.name in ok for g in circ.ops)


def strongest_witness(circ: Circuit):
    """Return (kind, detail) for the strongest sound whole-circuit witness
    available, or (None, reason) when none covers this circuit."""
    from .stabilizer import is_clifford

    if is_clifford(circ):
        return "stabilizer", {"note": "tableau equality, phase-exact up to "
                                      "global phase, any width"}
    if circ.num_qubits <= 8:
        return "unitary", {"note": "dense full-unitary re-computation"}
    if _is_phase_polynomial(circ):
        return "phase_polynomial", {
            "note": "parity/angle table equality, any width"}
    from .dd import DD_MAX_QUBITS
    if circ.num_qubits <= DD_MAX_QUBITS:
        return "dd", {"note": "decision-diagram overlap re-derived "
                              "independently by compactq-check (node-"
                              "budgeted; the checker declines, never "
                              "guesses)"}
    return None, ("circuit exceeds the dense prover and is neither "
                  "Clifford, phase-polynomial, nor DD-covered")


def build_certificate(original: Circuit, optimized: Circuit,
                      tol: float = 1e-7,
                      verification: Optional[dict] = None
                      ) -> tuple[Optional[dict], str]:
    """Produce the strongest whole-circuit certificate for the pair, or
    (None, reason) when no sound witness covers it.

    Soundness re-check is performed with compactq's provers here; the
    independent `compactq-check` program re-derives everything from the
    gate lists stored in the witness — no trusted numbers."""
    import compactq
    from .verify import verify
    from .stabilizer import is_clifford

    if original.num_qubits != optimized.num_qubits:
        return None, "qubit count mismatch"
    kind, _detail = strongest_witness(optimized)
    if kind is None:
        return None, "no sound whole-circuit witness for this width/family"

    v = verification if verification is not None else \
        verify(original, optimized)
    if v.get("equivalent") is not True:
        return None, f"verification did not confirm equivalence: {v}"

    # for stabilizer/phase-polynomial witnesses, confirm the family claim
    if kind == "stabilizer" and not (is_clifford(original)
                                     and is_clifford(optimized)):
        return None, "stabilizer witness requested but circuits are not Clifford"
    if kind == "phase_polynomial" and not _is_phase_polynomial(optimized):
        return None, "phase-polynomial witness requested but output is not diagonal+CNOT"

    witness = {"kind": kind, "n_qubits": original.num_qubits,
               "tolerance": tol}
    cert = {
        "cert_version": CERT_VERSION,
        "composition": "whole_circuit_v1",
        "producer": {"name": "compactq", "version": compactq.__version__},
        "input": {"qasm": to_qasm(original), "sha256": _qasm_sha(original)},
        "output": {"qasm": to_qasm(optimized), "sha256": _qasm_sha(optimized)},
        "n_qubits": original.num_qubits,
        "global_phase_convention": "equal_up_to_global_phase",
        "output_permutation": None,
        "witness": witness,
        "rewrite_trace": [],
        "claims": {"exact": True, "approximate": False, "permuted": False},
    }
    return cert, ""


def _piece_record(sub: Circuit) -> dict:
    return {"qasm": to_qasm(sub), "sha256": _qasm_sha(sub),
            "n_qubits": sub.num_qubits}


def build_certificate_blocks(original: Circuit, optimized: Circuit,
                             tol: float = 1e-7,
                             verification: Optional[dict] = None
                             ) -> tuple[Optional[dict], str]:
    """Disjoint-blocks compositional certificate
    (`composition: disjoint_blocks_v1`).

    Soundness (mirrors compactq.compositional): if the active qubits of
    both circuits partition into the SAME disjoint components, unitaries
    factorize over components, so per-block equivalence (up to global
    phase) implies whole-circuit equivalence.  Every block pair must sit
    within the dense prover's limit — no 2^total unitary is ever built.
    Returns (None, reason) when the structure does not apply, a piece is
    too wide, or any piece pair fails the in-product soundness re-check.
    """
    import compactq
    from .verify import verify
    from .compositional import active_components
    from .equivalence import check_equivalent, _MAX_QUBITS

    if original.num_qubits != optimized.num_qubits:
        return None, "qubit count mismatch"
    ca = active_components(original)
    cb = active_components(optimized)
    if len(ca) < 2 or len(cb) < 2:
        return None, "fewer than two disjoint components"
    if [qs for qs, _ in ca] != [qs for qs, _ in cb]:
        return None, "component structures differ between the circuits"
    if any(len(qs) > min(_MAX_QUBITS, _COMPOSITION_BLOCK_LIMIT)
           for qs, _ in ca):
        return None, "a component exceeds the per-piece dense limit"

    v = verification if verification is not None else \
        verify(original, optimized)
    if v.get("equivalent") is not True:
        return None, f"verification did not confirm equivalence: {v}"
    for (_qs, sub_a), (_qs2, sub_b) in zip(ca, cb):
        if not check_equivalent(sub_a, sub_b, tol):
            return None, "a block pair failed the in-product proof"

    cert = {
        "cert_version": CERT_VERSION,
        "composition": "disjoint_blocks_v1",
        "producer": {"name": "compactq", "version": compactq.__version__},
        "input": {"qasm": to_qasm(original), "sha256": _qasm_sha(original)},
        "output": {"qasm": to_qasm(optimized), "sha256": _qasm_sha(optimized)},
        "n_qubits": original.num_qubits,
        "global_phase_convention": "equal_up_to_global_phase",
        "output_permutation": None,
        "witness": {"kind": "composition",
                    "n_qubits": original.num_qubits,
                    "tolerance": tol},
        "blocks": {"original": [_piece_record(sub) for _qs, sub in ca],
                   "output": [_piece_record(sub) for _qs, sub in cb]},
        "rewrite_trace": [],
        "claims": {"exact": True, "approximate": False, "permuted": False},
    }
    return cert, ""


def build_certificate_segments(original: Circuit, optimized: Circuit,
                               cuts_a, cuts_b,
                               tol: float = 1e-7) -> tuple[Optional[dict], str]:
    """Sequential-segments certificate
    (`composition: sequential_segments_v1`).

    Soundness (mirrors compactq.compositional.verify_segmented): for any
    valid consecutive segmentations of the two circuits, per-segment
    equivalence (up to global phase) implies whole-circuit equivalence.
    Segments are dense-proven on their ACTIVE qubits only.  Cuts are the
    caller's rewrite correspondence (e.g. layer boundaries).  Returns
    (None, reason) when verify_segmented declines or a soundness
    re-check fails.
    """
    import compactq
    from .verify import verify
    from .compositional import verify_segmented

    if original.num_qubits != optimized.num_qubits:
        return None, "qubit count mismatch"
    v = verify(original, optimized)
    if v.get("equivalent") is not True:
        return None, f"verification did not confirm equivalence: {v}"
    seg = verify_segmented(original, optimized, cuts_a, cuts_b, tol)
    if seg is None:
        return None, ("segmentation not verifiable (width or structure "
                      "out of the segmented prover's scope)")
    if seg.get("equivalent") is not True:
        return None, "a segment pair failed the in-product proof"

    def slices(ops, cuts):
        bounds = [0] + sorted(int(x) for x in cuts) + [len(ops)]
        return [Circuit(original.num_qubits, ops[b:e])
                for b, e in zip(bounds, bounds[1:])]

    segs_a = slices(original.ops, cuts_a)
    segs_b = slices(optimized.ops, cuts_b)
    cert = {
        "cert_version": CERT_VERSION,
        "composition": "sequential_segments_v1",
        "producer": {"name": "compactq", "version": compactq.__version__},
        "input": {"qasm": to_qasm(original), "sha256": _qasm_sha(original)},
        "output": {"qasm": to_qasm(optimized), "sha256": _qasm_sha(optimized)},
        "n_qubits": original.num_qubits,
        "global_phase_convention": "equal_up_to_global_phase",
        "output_permutation": None,
        "witness": {"kind": "composition",
                    "n_qubits": original.num_qubits,
                    "tolerance": tol},
        "segments": {"original": [_piece_record(s) for s in segs_a],
                     "output": [_piece_record(s) for s in segs_b]},
        "rewrite_trace": [],
        "claims": {"exact": True, "approximate": False, "permuted": False},
    }
    return cert, ""


def optimize_with_certificate(circ: Circuit) -> "CertResult":
    """Optimize `circ` and return a CertResult carrying the optimized
    circuit plus its certificate (or the reason none was issued).

    Certificate ladder: whole-circuit witness first; if no whole-circuit
    witness covers the circuit, a disjoint-blocks compositional
    certificate (each block dense-proven); sequential-segment
    certificates are available explicitly via
    `build_certificate_segments` (they need segment correspondence
    which the optimizer does not track).

    Backward-compat note: `optimize_search(circ)` keeps its exact
    historical behavior and return type; this opt-in wrapper adds the
    certificate path."""
    from .search import optimize_search
    from .verify import verify

    opt = optimize_search(circ)
    v = verify(circ, opt)
    if v.get("equivalent") is not True:
        return CertResult(circuit=circ, certificate=None,
                          reason="verification declined — original returned")
    cert, reason = build_certificate(circ, opt, verification=v)
    if cert is not None:
        return CertResult(circuit=opt, certificate=cert, reason=reason)
    cert, reason = build_certificate_blocks(circ, opt, verification=v)
    if cert is not None:
        return CertResult(circuit=opt, certificate=cert, reason=reason)
    return CertResult(circuit=opt, certificate=None, reason=reason)
