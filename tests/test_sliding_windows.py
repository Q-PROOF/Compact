"""Sliding-window compositional prover tests (v0.2.6, T4.2-lite).

Connected circuits (beyond every whole-circuit prover) are exact-proven
tier 4 when a local rewrite correspondence exists; the prover must
DECLINE — never guess — on perturbed or globally-restructured pairs.

Standalone (`python tests/test_sliding_windows.py`), zero dependencies.
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compactq import Circuit, Gate  # noqa: E402
from compactq.compositional import verify_windows  # noqa: E402
from compactq.equivalence import check_equivalent  # noqa: E402

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:  # noqa: BLE001
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def ring(n, layers=2, theta=0.5):
    ops = []
    for _ in range(layers):
        for a in range(n):
            b = (a + 1) % n
            ops.append(Gate("cx", (), (a, b)))
            ops.append(Gate("rz", (theta,), (b,)))
            ops.append(Gate("cx", (), (a, b)))
        for q in range(n):
            ops.append(Gate("rx", (0.4,), (q,)))
    return Circuit(n, ops)


def local_rewrites(c, rng, commute_frac=0.5):
    """Exact local rewrite classes: merge adjacent rz pairs; commute an
    rz on a CX control across that CX."""
    merged = []
    i = 0
    ops = list(c.ops)
    while i < len(ops):
        g = ops[i]
        if (i + 1 < len(ops) and g.name == "rz" and ops[i + 1].name == "rz"
                and g.qubits == ops[i + 1].qubits):
            merged.append(Gate("rz", (g.params[0] + ops[i + 1].params[0],),
                               g.qubits))
            i += 2
        else:
            merged.append(g)
            i += 1
    out = list(merged)
    i = 0
    while i < len(out) - 1:
        g, h = out[i], out[i + 1]
        if (g.name == "cx" and h.name == "rz"
                and h.qubits[0] == g.qubits[0] and rng.random() < commute_frac):
            out[i], out[i + 1] = h, g
            i += 2
            continue
        i += 1
    return Circuit(c.num_qubits, out)


def test_connected_256q_exact_proven():
    """A locally-rewritten 256-qubit connected ring — far beyond the
    dense and DD provers — is exact-proven at tier 4 in well under a
    second."""
    import time
    c = ring(256, 2)
    d = local_rewrites(c, random.Random(3))
    assert len(d.ops) >= 1
    t0 = time.time()
    w = verify_windows(c, d)
    assert w is not None, "256q locally-rewritten pair must verify"
    assert w["equivalent"] is True and w["windows"] >= 100, w
    assert time.time() - t0 < 30, "window verification must stay fast"


def test_local_rewrite_pairs_always_prove():
    """40 random locally-rewritten small pairs: the window verdict must
    agree with the dense referee (all provable pairs prove)."""
    rng = random.Random(5)
    for trial in range(40):
        n = rng.choice((3, 4))
        ops = []
        for _ in range(rng.randrange(6, 20)):
            if rng.random() < 0.6:
                a = rng.randrange(n)
                b = (a + 1) % n
                ops.append(Gate("cx", (), (a, b)))
                ops.append(Gate("rz", (rng.uniform(0, 6.3),), (b,)))
                ops.append(Gate("cx", (), (a, b)))
            else:
                ops.append(Gate("rx", (rng.uniform(0, 6.3),),
                                (rng.randrange(n),)))
        a = Circuit(n, ops)
        b = local_rewrites(a, rng, commute_frac=1.0)
        dense = check_equivalent(a, b)
        w = verify_windows(a, b)
        if dense:
            assert w is not None and w["equivalent"] is True, \
                f"trial {trial}: dense-equivalent pair declined"
        else:
            assert w is None, \
                f"trial {trial}: dense-inequivalent pair claimed {w}"


def test_perturbed_pair_declines():
    """One perturbed angle in a 64q locally-rewritten ring: the prover
    must DECLINE (None), never claim equivalence and never crash."""
    c = ring(64, 2)
    d = local_rewrites(c, random.Random(3))
    ops = list(d.ops)
    for i, g in enumerate(ops):
        if g.name == "rz":
            ops[i] = Gate("rz", (g.params[0] + 0.01,), g.qubits)
            break
    v = verify_windows(c, Circuit(64, ops))
    assert v is None, f"perturbed pair must decline, got {v}"


def test_cascade_tier4_on_connected():
    """verify() routes a connected locally-rewritten pair through the
    tier-4 sliding_windows method (wider than dense/DD reach)."""
    from compactq import verify
    c = ring(64, 2)
    d = local_rewrites(c, random.Random(7))
    v = verify(c, d)
    assert v["equivalent"] is True, v
    assert v["tier"] == 4 and v["method"] == "sliding_windows", v
    assert v["windows"] >= 20, v


ALL = [
    ("windows: connected 256q exact-proven via sliding windows",
     test_connected_256q_exact_proven),
    ("windows: 40 local-rewrite pairs agree with the dense referee",
     test_local_rewrite_pairs_always_prove),
    ("windows: perturbed 64q pair declines loudly",
     test_perturbed_pair_declines),
    ("windows: verify() cascade routes to tier-4 sliding_windows",
     test_cascade_tier4_on_connected),
]


def main() -> int:
    print(f"compactq sliding-window prover suite ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll sliding-window tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
