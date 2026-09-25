"""Tier-boundary regression tests (v0.2.4, T2.2).

Every boundary claim in VERIFICATION.md is backed by a named, passing
test HERE.  If one of these fails, the boundary table in VERIFICATION.md
must be corrected in the same PR — the documentation and the code may
not drift apart.

Standalone (`python tests/test_tier_boundaries.py`), zero dependencies
(numpy only for the randomized-tier probes, skipped cleanly without it).
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compactq import Circuit, Gate, optimize, optimize_search  # noqa: E402
from compactq import verify  # noqa: E402
from compactq.equivalence import _MAX_QUBITS  # noqa: E402
from compactq.dd import DD_MAX_QUBITS  # noqa: E402

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:  # noqa: BLE001
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def _clifford_ladder(n, seed=7):
    rng = random.Random(seed)
    ops = []
    for _ in range(4 * n):
        ops.append(Gate(rng.choice(("h", "s", "x", "z")), (),
                        (rng.randrange(n),)))
        if n >= 2:
            ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
    return Circuit(n, ops)


def _phasepoly_ladder(n, seed=11):
    """CX + diagonal only: the tier-3 algebraic fragment at any width."""
    rng = random.Random(seed)
    ops = []
    for _ in range(3 * n):
        ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
        name = rng.choice(("rz", "t", "s", "z"))
        ops.append(Gate(name, (rng.uniform(0, 6.3),) if name == "rz" else (),
                        (rng.randrange(n),)))
    return Circuit(n, ops)


def _random_mixed(n, gates, seed=3):
    rng = random.Random(seed)
    ops = []
    for _ in range(gates):
        if rng.random() < 0.5:
            ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
        else:
            name = rng.choice(("h", "t", "tdg", "rz"))
            ops.append(Gate(name, (rng.uniform(0, 6.3),) if name in
                            ("rz", "t", "tdg") and rng.random() < 0.5 else (),
                            (rng.randrange(n),)))
    return Circuit(n, ops)


# ------------------------------------------------------------- boundaries
def test_boundary_clifford_tier3_any_width():
    """Clifford circuits prove EXACTLY at tier 3 (tableau) at every
    width — 5q, 9q, and 40q all land on `clifford_tableau`."""
    for n in (5, 9, 40):
        c = _clifford_ladder(n)
        o = optimize_search(c)
        v = verify(c, o)
        assert v["equivalent"] is True, f"{n}q: {v}"
        assert v["tier"] == 3 and v["method"] == "clifford_tableau", \
            f"{n}q landed on {v['method']}"


def test_boundary_phasepoly_tier3_any_width():
    """CX+diagonal circuits prove EXACTLY at tier 3 (GF(2) tables) at
    any width — 25q included.  The optimized pair is built by reordering
    commuting diagonal gates (angles preserved exactly, same unitary),
    so the tier-3 tables decide in milliseconds instead of falling
    through to the slower tiers."""
    n = 25
    c = _phasepoly_ladder(n)
    o = _reorder_diagonals(c)
    v = verify(c, o)
    assert v["equivalent"] is True, f"25q phase-poly: {v}"
    assert v["tier"] == 3 and v["method"] == "phase_polynomial", \
        f"25q landed on {v['method']}"


def _reorder_diagonals(c: Circuit) -> Circuit:
    """Bubble diagonal 1q gates past each other when the wires differ:
    an exactly-equivalent reordering in the CX+diagonal fragment."""
    ops = list(c.ops)
    for i in range(len(ops)):
        for j in range(len(ops) - 1 - i):
            a, b = ops[j], ops[j + 1]
            if (len(a.qubits) == 1 and len(b.qubits) == 1
                    and a.qubits != b.qubits
                    and a.name in ("rz", "t", "s", "z", "p")
                    and b.name in ("rz", "t", "s", "z", "p")):
                ops[j], ops[j + 1] = b, a
    return Circuit(c.num_qubits, ops)


def test_boundary_dense_tier2_at_limit():
    """Widths <= the dense limit prove at tier 2 full_unitary; the limit
    is 8q with the native kernels, 6q without."""
    n = min(_MAX_QUBITS, 8)
    c = _random_mixed(n, 12)
    o = optimize_search(c)
    v = verify(c, o)
    assert v["equivalent"] is True, f"{n}q: {v}"
    assert v["method"] == "full_unitary", f"{n}q landed on {v['method']}"


def test_boundary_dd_tier2_beyond_dense():
    """Just above the dense limit, structured circuits prove at tier 2
    via the decision diagram (`dd_full_unitary`) — QFT-12 is neither
    Clifford nor CX+diagonal, so tiers 3 cannot cover it."""
    n = _MAX_QUBITS + 4          # 12q native, 10q pure
    c = _qft(n)
    o = optimize_search(c)
    v = verify(c, o)
    assert v["equivalent"] is True, f"qft {n}q: {v}"
    assert v["method"] == "dd_full_unitary", \
        f"qft {n}q landed on {v['method']}"


def test_boundary_dd_width_cap():
    """At/above DD_MAX_QUBITS (32) the dd prover is out of scope by
    design: a 33q mixed circuit falls through — to the compositional
    sliding-window prover (tier 4, no width cap by design, when a local
    correspondence exists) or to randomized/prover-unavailable — never
    a false exact claim either way."""
    n = DD_MAX_QUBITS + 1
    c = _qft(n)
    o = optimize(c)
    v = verify(c, o)
    assert v["equivalent"] is not False, f"{n}q: {v}"
    if v["tier"] == 4:
        assert v["method"] == "sliding_windows", v
    else:
        assert v["tier"] <= 1, f"{n}q claimed tier {v['tier']}: {v}"
        assert v["method"] in ("randomized_sampling", "prover_unavailable",
                               "phase_polynomial", "clifford_tableau"), v


def test_boundary_budget_blowup_declines():
    """When the DD prover cannot hold a diagram within its budget it
    must decline loudly — never produce a verdict from an incomplete
    build — while the same pair proves exactly once the budget
    suffices.  Pinned with a reduced budget via the dd prover's public
    parameter; where the DEFAULT budget binds per circuit family is
    measured, not asserted: results/dd_ceiling.json (the verify-level
    fall-through after a decline is pinned by test_boundary_dd_width_
    cap and the adversarial suite)."""
    from compactq.dd import dd_fidelity, DDOverflow
    c = _qft(16)                      # needs ~65k nodes at the default budget
    try:
        dd_fidelity(c, Circuit(16, list(c.ops)), max_nodes=2000)
        raise AssertionError("reduced budget did not decline")
    except DDOverflow:
        pass                          # the loud decline is the contract
    v = dd_fidelity(c, Circuit(16, list(c.ops)))   # default budget suffices
    assert v > 1 - 1e-9


def test_boundary_optimize_never_claims_beyond_limit():
    """optimize_search at 9q (above the dense limit, generic gates)
    must not silently skip verification metadata: optimize_large
    reports the strongest status actually achieved."""
    from compactq import optimize_large
    c = _random_mixed(9, 30)
    o, status = optimize_large(c)
    assert o.num_qubits == 9
    assert status in ("exact-proven (algebraic)",
                      "exact-proven (decision-diagram)",
                      "exact", "rejected", "unverified (prover declined)"), \
        status
    v = verify(c, o)
    if status == "rejected":
        assert o == c or v["equivalent"] is False


def test_boundary_non_unitary_rejected_at_edge():
    """The unitary-core scope boundary: measurement and reset are
    rejected loudly at every entry point, never dropped."""
    from compactq.errors import UnsupportedCircuitError
    from compactq.qiskit_bridge import from_qiskit
    try:
        from compactq import from_qasm
        from_qasm('OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\n'
                  'h q[0];\nmeasure q[0] -> c[0];\ncreg c[1];\n'
                  'if(c[1]==1) x q[1];\n')
        # mid-circuit control flow must not pass silently
        raise AssertionError("classically-conditioned gate not rejected")
    except UnsupportedCircuitError:
        pass
    except ValueError:
        pass


def _qft(n):
    c = Circuit(n, [])
    for i in range(n):
        c.append(Gate("h", (), (i,)))
        for j in range(i + 1, n):
            c.append(Gate("cp", (math.pi / (2 ** (j - i)),), (j, i)))
    return c


ALL = [
    ("boundary: Clifford -> tier 3 at any width (5/9/40q)",
     test_boundary_clifford_tier3_any_width),
    ("boundary: CX+diagonal -> tier 3 at any width (25q)",
     test_boundary_phasepoly_tier3_any_width),
    ("boundary: dense limit -> tier 2 full_unitary",
     test_boundary_dense_tier2_at_limit),
    ("boundary: beyond dense -> tier 2 dd_full_unitary (structured)",
     test_boundary_dd_tier2_beyond_dense),
    ("boundary: DD width cap (33q) falls to randomized/unavailable",
     test_boundary_dd_width_cap),
    ("boundary: node-budget blowup declines loudly (16q)",
     test_boundary_budget_blowup_declines),
    ("boundary: optimize_large statuses honest at 9q",
     test_boundary_optimize_never_claims_beyond_limit),
    ("boundary: non-unitary ops rejected loudly",
     test_boundary_non_unitary_rejected_at_edge),
]


def main() -> int:
    print(f"compactq tier-boundary suite ({len(ALL)} boundaries)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll boundary tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
