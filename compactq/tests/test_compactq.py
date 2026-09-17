"""pytest-compatible tests for compactq."""
import pytest
from compactq import Circuit, Gate, optimize, optimize_search, from_qasm, to_qasm, benchmarks
from compactq.equivalence import check_equivalent


def test_u3_import():
    c = from_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[1];\nu3(0.7,1.3,-0.4) q[0];")
    assert len(c.ops) == 3  # rz, ry, rz


def test_ccx_import():
    c = from_qasm("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nccx q[0],q[1],q[2];")
    assert len(c.ops) == 15  # exact CCX decomposition


def test_optimize_preserves_2q_count():
    c = benchmarks.qft(4)
    o = optimize(c)
    assert o.two_qubit_count() <= c.two_qubit_count()


@pytest.mark.parametrize("name", list(benchmarks.default_suite().keys()))
def test_suite_verified(name):
    c = benchmarks.default_suite()[name]
    o = optimize(c)
    assert check_equivalent(c, o)


def test_optimize_search_no_regression():
    for name, c in benchmarks.default_suite().items():
        o = optimize_search(c)
        assert check_equivalent(c, o)
        assert len(o.ops) <= len(c.ops)
