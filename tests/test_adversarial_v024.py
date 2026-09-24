"""Adversarial soundness probes, v0.2.4 red-team (T1.3).

Dedicated pre-release red-team of the prover stack: every probe must
return a correct verdict or decline loudly.  Finds the class of bug a
random-circuit property sweep misses — degenerate parameters, global
phase traps, tier-straddling families, budget boundaries, malformed
input.  Findings fixed this release: same-wire two-qubit gates were
silently accepted by the IR; qelib1's `ch` gate failed QASM import.

Standalone (`python tests/test_adversarial_v024.py`), zero dependencies.
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compactq import Circuit, Gate, optimize, optimize_deep, optimize_search  # noqa: E402
from compactq.equivalence import check_equivalent  # noqa: E402

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:  # noqa: BLE001
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def _large(c):
    from compactq import optimize_large
    return optimize_large(c)


# ------------------------------------------------------------------ tests
def test_same_wire_gates_rejected():
    """cx(0,0) has no well-defined unitary; the IR must reject it at
    construction (found via adversarial testing, v0.2.4)."""
    for name in ("cx", "cz", "swap", "cp"):
        kw = ((0.3,) if name == "cp" else ())
        try:
            Gate(name, kw, (0, 0))
            raise AssertionError(f"{name} on (0,0) was accepted")
        except ValueError:
            pass
    try:
        Gate("mcx", (), (1, 1))
        raise AssertionError("mcx on duplicate wires accepted")
    except ValueError:
        pass


def test_ch_qasm_import_exact():
    """qelib1 `ch` imports to the exact controlled-H (adversarial
    testing find, v0.2.4): dense unitary must match to 1e-13, both wire
    orders."""
    from compactq import from_qasm
    from compactq.equivalence import unitary
    s2 = 1 / math.sqrt(2)
    for ctrl in (0, 1):
        tgt = 1 - ctrl
        c = from_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"
                      f"qreg q[2];\nch q[{ctrl}], q[{tgt}];\n")
        U = [x for row in unitary(c) for x in row]
        d = 4
        ref = [0.0j] * (d * d)
        for col in range(d):              # column = input basis state
            cbit = (col >> ctrl) & 1
            tbit = (col >> tgt) & 1
            if cbit == 0:
                ref[col * d + col] += 1.0
            else:
                for new_t, amp in ((0, s2), (1, -s2 if tbit else s2)):
                    row = (col & ~(1 << tgt)) | (new_t << tgt)
                    ref[row * d + col] += amp
        tr = sum(u.conjugate() * r for u, r in zip(U, ref))
        fid = abs(tr) / d
        assert fid > 1 - 1e-12, f"ch fidelity {fid} (ctrl={ctrl})"


def test_global_phase_traps():
    """y.y = -I must collapse to empty; x/y orderings equal up to phase."""
    yy = Circuit(1, [Gate("y", (), (0,)), Gate("y", (), (0,))])
    o = optimize(yy)
    assert len(o) == 0, f"y.y -> {o.stats()}"
    a = Circuit(2, [Gate("x", (), (0,)), Gate("y", (), (1,))])
    b = Circuit(2, [Gate("y", (), (1,)), Gate("x", (), (0,))])
    assert check_equivalent(a, b)
    o = optimize_search(a)
    assert check_equivalent(b, o)


def test_identity_chains_collapse():
    """Gate + exact inverse, interleaved: must optimize to the empty
    circuit, exactly, on 60 random chains."""
    rng = random.Random(2024)
    for trial in range(60):
        n = rng.choice((1, 2, 3))
        ops = []
        for _ in range(rng.randrange(4, 30)):
            g = rng.choice(("h", "x", "y", "z", "s", "sdg", "t", "tdg"))
            ops.append(Gate(g, (), (rng.randrange(n),)))
        inv = []
        for g in reversed(ops):
            if g.name in ("s", "t"):
                inv.append(Gate(g.name + "dg", (), g.qubits))
            elif g.name in ("sdg", "tdg"):
                inv.append(Gate(g.name[:-2], (), g.qubits))
            else:
                inv.append(Gate(g.name, (), g.qubits))
        c = Circuit(n, ops + inv)
        o = optimize(c)
        assert len(o) == 0, f"trial {trial}: {o.stats()}"


def test_degenerate_angles():
    """0, -0, 2pi, -2pi, 1e-18, 100pi: optimize must preserve the
    unitary exactly and never grow."""
    for ang in (0.0, -0.0, 1e-18, -1e-18, 2 * math.pi, 2 * math.pi + 1e-18,
                -2 * math.pi, 4 * math.pi, 100 * math.pi, math.pi, -math.pi):
        c = Circuit(2, [Gate("rz", (ang,), (0,)),
                        Gate("p", (ang,), (1,)),
                        Gate("rx", (ang,), (0,))])
        o = optimize(c)
        assert check_equivalent(c, o), f"angle {ang}"
        assert len(o) <= len(c)


def test_tier_straddling_clifford_t():
    """9q Clifford+T with H gates: not Clifford, not CX+diagonal, above
    the dense ceiling — must prove at the dd tier (or decline), never
    claim an algebraic tier that cannot cover it."""
    from compactq import verify
    from compactq.cert import optimize_with_certificate
    rng = random.Random(5)
    ops = []
    for _ in range(24):
        ops.append(Gate(rng.choice(("h", "t", "tdg", "s")), (),
                        (rng.randrange(9),)))
        ops.append(Gate("cx", (), tuple(rng.sample(range(9), 2))))
    c = Circuit(9, ops)
    o, status = _large(c)
    v = verify(c, o)
    assert v["equivalent"] is True, f"9q H+T unproven: {v} ({status})"
    assert v["method"] == "dd_full_unitary", \
        f"expected dd tier, got {v['method']}"
    res = optimize_with_certificate(c)
    assert res.certificate is not None, res.reason
    assert res.certificate["witness"]["kind"] == "dd"


def test_dd_budget_boundary_declines_loudly():
    """Wide random circuits stress the DD budget boundary: whatever the
    outcome, the cascade must be sound — exact tiers only when a prover
    actually covered the circuit (post-v0.2.4 max-weight normalization
    some 16q random circuits now prove at the dd tier), randomized tier
    when it declines, and never a wrong verdict or a crash."""
    from compactq import verify
    rng = random.Random(77)
    for trial in range(1):
        ops = []
        for _ in range(60):
            ops.append(Gate("cx", (), tuple(rng.sample(range(16), 2))))
            ops.append(Gate("rz", (rng.uniform(0, 6.3),), (rng.randrange(16),)))
        c = Circuit(16, ops)
        o, status = _large(c)
        assert o.num_qubits == 16
        v = verify(c, o)
        assert v["equivalent"] is not False, f"trial {trial}: {v}"
        if v["equivalent"] is True:
            assert v["tier"] >= 1, f"trial {trial}: {v}"
        assert "unverified" not in status or status == "unverified"


def test_cross_prover_verdict_agreement():
    """Dense, decision-diagram and phase-polynomial provers must agree on
    clearly-equal, clearly-different and noise-level-perturbed pairs."""
    from compactq.dd import check_equivalent_dd
    from compactq.phasepoly import phasepoly_equal
    rng = random.Random(31)
    for trial in range(120):
        n = rng.choice((3, 4))
        ops = []
        for _ in range(rng.randrange(6, 24)):
            ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
            ops.append(Gate("rz", (rng.uniform(0, 6.3),), (rng.randrange(n),)))
        c1 = Circuit(n, ops)
        mode = rng.random()
        if mode < 0.5:
            c2 = Circuit(n, list(ops))
        elif mode < 0.8:  # clearly different (1e-3 >> every tolerance)
            k = rng.randrange(len(ops))
            g = ops[k]
            newp = (g.params[0] + 1e-3,) if g.params else ()
            ops2 = list(ops)
            ops2[k] = Gate(g.name, newp, g.qubits)
            c2 = Circuit(n, ops2)
        else:            # 1e-11 noise: same under every prover tolerance
            k = rng.randrange(len(ops))
            g = ops[k]
            newp = (g.params[0] + 1e-11,) if g.params else ()
            ops2 = list(ops)
            ops2[k] = Gate(g.name, newp, g.qubits)
            c2 = Circuit(n, ops2)
        dense = check_equivalent(c1, c2)
        dd = check_equivalent_dd(c1, c2)
        assert dense == dd, f"trial {trial}: dense={dense} dd={dd}"
        pp = phasepoly_equal(c1, c2)
        if pp is not None:
            assert pp == dense, f"trial {trial}: dense={dense} pp={pp}"


def test_verify_argument_symmetry():
    """verify(a, b) and verify(b, a) must reach the same verdict."""
    from compactq import verify
    rng = random.Random(9)
    for _ in range(10):
        n = rng.choice((2, 3))
        ops = []
        for _ in range(rng.randrange(3, 15)):
            if rng.random() < 0.5:
                ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
            else:
                ops.append(Gate(rng.choice(("h", "t", "rz")),
                                (rng.uniform(0, 6.3),) if rng.random() < 0.5
                                else (), (rng.randrange(n),)))
        c = Circuit(n, ops)
        o = optimize_search(c)
        va = verify(c, o)
        vb = verify(o, c)
        assert va["equivalent"] == vb["equivalent"], f"{va} vs {vb}"


def test_never_grow_on_nasty_structures():
    """CP ladders and 200-gate 1q runs: no objective order may grow."""
    rng = random.Random(4242)
    families = []
    for n in (2, 3, 4):
        ops = []
        for j in range(20):
            ops.append(Gate("cp", (rng.uniform(0, 6.3),), (j % n, (j + 1) % n)))
            ops.append(Gate("rz", (rng.uniform(0, 6.3),), (j % n,)))
        families.append(Circuit(n, ops))
    for _ in range(3):
        ops = []
        for _ in range(200):
            ops.append(Gate(rng.choice(("h", "x", "rz", "t", "sdg", "sx")),
                            (rng.uniform(0, 6.3),) if rng.random() < 0.3
                            else (), (0,)))
        families.append(Circuit(1, ops))
    for c in families:
        for obj in ("2q", "depth", "gate_count"):
            outs = [optimize(c, objective=obj),
                    optimize_deep(c, objective=obj), optimize_search(c)]
            for o in outs:
                got = (o.two_qubit_count() if obj == "2q"
                       else o.depth() if obj == "depth" else len(o))
                base = (c.two_qubit_count() if obj == "2q"
                        else c.depth() if obj == "depth" else len(c))
                assert got <= base, f"{obj} grew: {base} -> {got}"


def test_certificate_edges():
    """Empty circuit and phase-only collapse still certify soundly."""
    from compactq.cert import optimize_with_certificate
    res = optimize_with_certificate(Circuit(2, []))
    assert res.certificate is not None, res.reason
    yy = Circuit(3, [Gate("y", (), (0,)), Gate("y", (), (0,))])
    res2 = optimize_with_certificate(yy)
    assert res2.certificate is not None, res2.reason


def test_dd_self_comparison_mcx():
    """A circuit compared against its own copy must get self-fidelity 1
    from the dd prover at every width.  Found via v0.2.4 red-teaming:
    first-nonzero normalization amplified near-cancelling blocks and a
    tier-2 INEQUIVALENT verdict was returned for identical 8q MCX
    circuits.  Max-weight normalization fixes the arithmetic; the
    self-consistency norm guard turns any residue into a decline."""
    from compactq.dd import dd_fidelity, DDOverflow
    for k in range(2, 9):
        c = Circuit(k + 1, [Gate("mcx", (), tuple(range(k + 1)))])
        try:
            v = dd_fidelity(c, Circuit(k + 1, list(c.ops)))
            assert v > 1 - 1e-9, f"mcx {k}c: self-fidelity {v}"
        except DDOverflow:
            pass  # loud decline is acceptable, a wrong verdict is not
    for n in (6, 8, 10):
        ops = [Gate("h", (), (q,)) for q in range(n)]
        ops.append(Gate("mcx", (), tuple(range(n))))
        c = Circuit(n, ops)
        v = dd_fidelity(c, Circuit(n, list(c.ops)))
        assert v > 1 - 1e-9, f"grover {n}q: self-fidelity {v}"


ALL = [
    ("adversarial: same-wire gates rejected at construction",
     test_same_wire_gates_rejected),
    ("adversarial: dd self-comparison on mcx chains (regression)",
     test_dd_self_comparison_mcx),
    ("adversarial: qelib1 ch import exact (both wire orders)",
     test_ch_qasm_import_exact),
    ("adversarial: global-phase traps (y.y collapse, x/y order)",
     test_global_phase_traps),
    ("adversarial: identity chains collapse exactly (60)",
     test_identity_chains_collapse),
    ("adversarial: degenerate angles preserved exactly",
     test_degenerate_angles),
    ("adversarial: 9q Clifford+T straddles to the dd tier",
     test_tier_straddling_clifford_t),
    ("adversarial: DD budget boundary declines loudly (16q)",
     test_dd_budget_boundary_declines_loudly),
    ("adversarial: cross-prover verdict agreement (120 pairs)",
     test_cross_prover_verdict_agreement),
    ("adversarial: verify argument symmetry",
     test_verify_argument_symmetry),
    ("adversarial: never-grow on nasty structures",
     test_never_grow_on_nasty_structures),
    ("adversarial: certificate edge cases", test_certificate_edges),
]


def main() -> int:
    print(f"compactq adversarial red-team suite ({len(ALL)} probes)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll adversarial probes passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
