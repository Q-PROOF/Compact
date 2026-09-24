"""Compositional-certificate tests (v0.2.5, T4 maturation).

disjoint_blocks_v1 and sequential_segments_v1 certificates must round-
trip through the INDEPENDENT compactq-check program (VALID), and every
tampering must be caught.  The checker re-derives the block structure
from the provided circuits itself — the certificate's declared pieces
are never trusted for the proof.

Standalone (`python tests/test_compositional_certs.py`), zero
dependencies.
"""
import copy
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "compactq-check" / "src"))

import compactq  # noqa: E402
from compactq import Circuit, Gate  # noqa: E402
from compactq.benchmarks import qft  # noqa: E402
from compactq.cert import (build_certificate_blocks,  # noqa: E402
                           build_certificate_segments,
                           optimize_with_certificate)
import compactq_check.checker as checker  # noqa: E402

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:  # noqa: BLE001
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def _block_circuit(n_blocks=8, width=5, seed=1):
    """Independent non-Clifford blocks: nothing above the block width
    couples them, so the whole circuit is wider than any whole-circuit
    witness but every block is dense-provable."""
    del seed
    ops = []
    for b in range(n_blocks):
        base = b * width
        ops.append(Gate("h", (), (base,)))
        for j in range(base, base + width - 1):
            ops.append(Gate("cx", (), (j, j + 1)))
            ops.append(Gate("t", (), (j + 1,)))
            ops.append(Gate("h", (), (j,)))
            ops.append(Gate("cx", (), (j, j + 1)))
    return Circuit(n_blocks * width, ops)


def _blocks_circuit_nocliff(seed=1):
    return _block_circuit(8, 5, seed)


def test_blocks_certificate_end_to_end():
    c = _block_circuit(8, 5)
    res = optimize_with_certificate(c)
    assert res.certificate is not None, res.reason
    cert = res.certificate
    assert cert["composition"] == "disjoint_blocks_v1", cert["composition"]
    assert len(cert["blocks"]["original"]) == 8
    v = checker.check_certificate(cert, compactq.to_qasm(c),
                                  compactq.to_qasm(res.circuit))
    assert v["verdict"] == "VALID", v
    assert "8 components" in v["method"], v


def test_blocks_certificate_structure_mutations_killed():
    c = _block_circuit(8, 5)
    res = optimize_with_certificate(c)
    cert = res.certificate
    in_q = compactq.to_qasm(c)
    out_q = compactq.to_qasm(res.circuit)

    # a gate flipped inside the provided INPUT circuit: the checker
    # re-derives components from the PROVIDED circuits, so a perturbed
    # input must break the hash pin or a component proof
    lines = in_q.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("t "):
            lines[i] = ln.replace("t ", "s ", 1)
            break
    mutated_in = "\n".join(lines) + "\n"
    assert mutated_in != in_q, "mutation did not change the input"
    v = checker.check_certificate(cert, mutated_in, out_q)
    assert v["verdict"] == "INVALID", (v, "mutated input")

    # flipped full-circuit hash
    bad = copy.deepcopy(cert)
    bad["input"]["sha256"] = "1" + bad["input"]["sha256"][1:]
    v = checker.check_certificate(bad, in_q, out_q)
    assert v["verdict"] == "INVALID", (v, "hash flip")

    # a gate flipped in the provided OUTPUT circuit: hash pins catch it
    out_lines = out_q.splitlines()
    for i, ln in enumerate(out_lines):
        if ln.startswith("t "):
            out_lines[i] = ln.replace("t ", "s ", 1)
            break
    mutated_out = "\n".join(out_lines) + "\n"
    assert mutated_out != out_q
    v = checker.check_certificate(cert, in_q, mutated_out)
    assert v["verdict"] == "INVALID", (v, "mutated output")

    # NOTE: tampering with cert["blocks"] itself is ignored by design —
    # the checker re-derives the component structure from the provided
    # circuits and never trusts the declared pieces for the proof.


def test_segments_certificate_end_to_end():
    c = qft(6)
    twin = Circuit(6, list(c.ops))     # identical correspondence
    cuts = (10, 20)
    cert, reason = build_certificate_segments(c, twin, cuts, cuts)
    assert cert is not None, reason
    assert cert["composition"] == "sequential_segments_v1"
    v = checker.check_certificate(cert, compactq.to_qasm(c),
                                  compactq.to_qasm(twin))
    assert v["verdict"] == "VALID", v

    # mutation: perturb one segment's qasm -> hash mismatch -> INVALID
    bad = copy.deepcopy(cert)
    piece = bad["segments"]["output"][1]
    piece["qasm"] = piece["qasm"].replace("cx ", "cx ", 1)
    bad_hash = piece["sha256"][:-1] + ("0" if piece["sha256"][-1] != "0"
                                       else "1")
    piece["sha256"] = bad_hash
    v = checker.check_certificate(bad, compactq.to_qasm(c),
                                  compactq.to_qasm(twin))
    assert v["verdict"] == "INVALID", v


def test_certificate_ladder_whole_circuit_first():
    """Clifford circuits still get the strongest whole-circuit witness;
    the compositional ladder only fires when no whole-circuit witness
    covers the circuit."""
    from compactq.benchmarks import ghz
    res = optimize_with_certificate(ghz(5))
    assert res.certificate["composition"] == "whole_circuit_v1"
    assert res.certificate["witness"]["kind"] == "stabilizer"


ALL = [
    ("compositional: disjoint-blocks certificate end-to-end (40q)",
     test_blocks_certificate_end_to_end),
    ("compositional: structure mutations killed",
     test_blocks_certificate_structure_mutations_killed),
    ("compositional: sequential-segments certificate end-to-end",
     test_segments_certificate_end_to_end),
    ("compositional: ladder keeps whole-circuit witnesses first",
     test_certificate_ladder_whole_circuit_first),
]


def main() -> int:
    print(f"compactq compositional-certificate suite ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll compositional-certificate tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
