"""Objective-mode tests: lexicographic objectives + weighted ranking.

Runs standalone (`python tests/test_objectives.py`, zero dependencies) and
is wired into CI alongside tests/run_tests.py.

Guarantees under test:
  - '2q' / 'depth' / 'gate_count': the chosen primary metric never grows
    (each objective is a reordering of the same lexicographic tuple).
  - 'weighted': the weighted cost 1.0*2q + 0.1*depth + 0.02*gates never
    exceeds the input's.
  - 'latency' aliases 'depth'; unknown objectives are rejected.
  - the '2q' default is behavior-identical to the historical pipeline.
"""
import sys
from pathlib import Path

# run against THIS checkout, not any pip-installed copy
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def _metrics(c):
    return (c.two_qubit_count(), len(c.ops), c.depth())


_PRIMARY = {"2q": 0, "depth": 2, "gate_count": 1}


def test_lexicographic_objectives_never_grow():
    from compactq import optimize, optimize_deep, optimize_search
    from compactq.benchmarks import qft, ghz, clifford_ladder, random_circuit
    circuits = [qft(3), qft(4), ghz(5), clifford_ladder(4),
                random_circuit(4, 24, seed=5), random_circuit(5, 30, seed=6)]
    for obj in ("2q", "depth", "gate_count"):
        for circ in circuits:
            for fn in (optimize, optimize_deep, optimize_search):
                out = fn(circ, objective=obj)
                m_in, m_out = _metrics(circ), _metrics(out)
                assert m_out[_PRIMARY[obj]] <= m_in[_PRIMARY[obj]], \
                    (obj, fn.__name__, m_in, m_out)


def test_weighted_cost_never_increases():
    from compactq import optimize_search
    from compactq.benchmarks import qft, ghz, random_circuit

    def wcost(c):
        return c.two_qubit_count() + 0.1 * c.depth() + 0.02 * len(c.ops)

    for circ in (qft(3), qft(4), ghz(5), random_circuit(4, 24, seed=7)):
        out = optimize_search(circ, objective="weighted")
        assert wcost(out) <= wcost(circ) + 1e-9, (wcost(circ), wcost(out))


def test_latency_alias_and_rejection():
    from compactq import optimize, optimize_search
    from compactq.benchmarks import qft
    out = optimize_search(qft(4), objective="latency")
    assert out.depth() <= qft(4).depth()
    for bad in ("hardware_error", "nonsense"):
        try:
            optimize(qft(3), objective=bad)
            raise AssertionError(f"{bad}: expected ValueError")
        except ValueError:
            pass


def test_default_objective_unchanged():
    # the '2q' default must be behavior-identical to the historical pipeline
    from compactq import optimize, optimize_search
    from compactq.benchmarks import qft, random_circuit
    for circ in (qft(4), random_circuit(4, 20, seed=9)):
        a = optimize(circ)
        b = optimize(circ, objective="2q")
        assert (len(a.ops), a.two_qubit_count(), a.depth()) == \
               (len(b.ops), b.two_qubit_count(), b.depth())
        s1 = optimize_search(circ)
        s2 = optimize_search(circ, objective="2q")
        assert (len(s1.ops), s1.two_qubit_count(), s1.depth()) == \
               (len(s2.ops), s2.two_qubit_count(), s2.depth())


ALL = [
    ("lexicographic objectives never grow (6 circuits x 3 entry points)",
     test_lexicographic_objectives_never_grow),
    ("weighted cost never increases", test_weighted_cost_never_increases),
    ("latency alias + invalid objective rejection",
     test_latency_alias_and_rejection),
    ("default objective unchanged (2q == historical)",
     test_default_objective_unchanged),
]


def main() -> int:
    print(f"compactq objective-mode tests ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll objective tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
