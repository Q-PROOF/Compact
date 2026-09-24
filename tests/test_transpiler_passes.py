"""Integration tests: compactq as passes inside existing toolchains
(v0.2.4, T3.1 + T3.2).

Acceptance contract (release plan):
  * `PassManager([CompactqPass()]).run(qc)` output matches the direct
    compactq call on the same input (gate-exact + unitary fidelity),
    and the verification verdict survives as inspectable pass output
    (property_set + circuit metadata), never silently dropped.
  * `CompactTketPass().apply(tk)` behaves identically through pytket,
    with `last_proof` carrying the verdict.
  * the pytket/qiskit bridges are exact in both directions.

Skips gracefully (exit 0) when qiskit/pytket are not installed — the
core package stays zero-dependency.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:  # noqa: BLE001
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


try:
    import numpy as np
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import Operator
    from qiskit.transpiler import PassManager
    HAVE_QISKIT = True
except Exception:
    HAVE_QISKIT = False

try:
    import pytket  # noqa: F401
    HAVE_TKET = True
except Exception:
    HAVE_TKET = False


# ------------------------------------------------------------------ qiskit
def test_qiskit_pass_matches_direct_call():
    if not HAVE_QISKIT:
        print("    (skipped: qiskit not installed)")
        return
    from compactq.plugins.qiskit_plugin import CompactqPass
    from compactq.qiskit_bridge import compactq_pass

    rng_cases = [_qiskit_case(n) for n in (3, 4)]
    for qc in rng_cases:
        via_pm = PassManager([CompactqPass()]).run(qc)
        direct = compactq_pass(qc)
        assert dict(via_pm.count_ops()) == dict(direct.count_ops()), \
            f"pass vs direct gate counts differ: {via_pm.count_ops()} " \
            f"vs {direct.count_ops()}"
        f = _fidelity(via_pm, direct)
        assert f > 1 - 1e-9, f"pass vs direct fidelity {f}"


def test_qiskit_pass_proof_metadata():
    if not HAVE_QISKIT:
        print("    (skipped: qiskit not installed)")
        return
    from compactq.plugins.qiskit_plugin import CompactqPass

    qc = _qiskit_case(4)
    pm = PassManager([CompactqPass()])
    out = pm.run(qc)
    ps = pm.property_set
    proof = ps.get("compactq_proof") if hasattr(pm, "property_set") else None
    assert proof is not None, "property_set['compactq_proof'] missing"
    assert proof["equivalent"] is True and proof["tier"] >= 2, proof
    meta = (out.metadata or {}).get("compactq")
    assert meta is not None, "output.metadata['compactq'] missing"
    assert meta["equivalent"] is True


def test_qiskit_pass_level0_and_approximate_paths():
    if not HAVE_QISKIT:
        print("    (skipped: qiskit not installed)")
        return
    from compactq.plugins.qiskit_plugin import CompactqPass

    qc = _qiskit_case(3)
    out0 = PassManager([CompactqPass(optimization_level=0)]).run(qc)
    f = _fidelity(out0, qc)
    assert f > 1 - 1e-9, f"level-0 pass fidelity {f}"
    outa = PassManager([CompactqPass(approximate=True,
                                     min_fidelity=0.999)]).run(qc)
    fa = _fidelity(outa, qc)
    assert fa >= 0.99, f"approximate pass fidelity {fa}"


def _qiskit_case(n):
    qc = QuantumCircuit(n)
    for q in range(n):
        qc.h(q)
    for q in range(n - 1):
        qc.cx(q, q + 1)
        qc.t(q)
        qc.cx(q, q + 1)
    qc.sdg(0)
    qc.cx(n - 1, 0)
    return qc


def _fidelity(a, b):
    na, nb = a.num_qubits, b.num_qubits
    assert na == nb
    va = Operator(a).data
    vb = Operator(b).data
    return float(abs(np.sum(np.conj(va) * vb)) / (1 << na))


# ------------------------------------------------------------------- pytket
def test_tket_pass_matches_direct_call():
    if not (HAVE_TKET and HAVE_QISKIT):
        print("    (skipped: pytket/qiskit not installed)")
        return
    from compactq.plugins.pytket_plugin import CompactTketPass
    from compactq.pytket_bridge import to_pytket, from_pytket
    from compactq import optimize_search
    from compactq.qiskit_bridge import to_qiskit
    from pytket.extensions.qiskit import tk_to_qiskit, qiskit_to_tk

    qc = _qiskit_case(4)
    tk = qiskit_to_tk(qc)
    out_tk = CompactTketPass().apply(tk)
    direct = optimize_search(__import__("compactq.qiskit_bridge",
                                        fromlist=["from_qiskit"])
                             .from_qiskit(qc))
    f = _fidelity(tk_to_qiskit(out_tk), to_qiskit(direct))
    assert f > 1 - 1e-9, f"tket pass vs direct fidelity {f}"


def test_tket_pass_proof_metadata_and_sequence():
    if not (HAVE_TKET and HAVE_QISKIT):
        print("    (skipped: pytket/qiskit not installed)")
        return
    from compactq.plugins.pytket_plugin import CompactTketPass
    from pytket.passes import SequencePass
    from pytket.extensions.qiskit import qiskit_to_tk

    tk = qiskit_to_tk(_qiskit_case(3))
    p = CompactTketPass()
    out = p.apply(tk)
    assert p.last_proof is not None, "last_proof missing"
    assert p.last_proof["equivalent"] is True, p.last_proof
    seq = SequencePass([p.as_tket_pass()])
    ok = seq.apply(tk)
    assert ok, "SequencePass application failed"
    assert p.last_proof["equivalent"] is True, "proof lost in SequencePass"


# ------------------------------------------------------------------ bridges
def test_pytket_bridge_exact_both_directions():
    if not (HAVE_TKET and HAVE_QISKIT):
        print("    (skipped: pytket/qiskit not installed)")
        return
    import math
    import random
    from compactq import Circuit, Gate
    from compactq.pytket_bridge import to_pytket, from_pytket
    from compactq.qiskit_bridge import to_qiskit
    from pytket.extensions.qiskit import tk_to_qiskit
    from pytket.circuit import Circuit as TC

    def f2(qc_a, qc_b):
        va, vb = Operator(qc_a).data, Operator(qc_b).data
        return float(abs(np.sum(np.conj(va) * vb)) / va.shape[0])

    # tket-native constructs -> compactq
    t = TC(3)
    t.TK1(0.3, 0.4, 0.5, 0)
    t.CCX(0, 1, 2)
    t.CSWAP(2, 1, 0)
    t.CRz(0.7, 0, 1)
    t.SWAP(0, 2)
    f = f2(tk_to_qiskit(t), to_qiskit(from_pytket(t)))
    assert f > 1 - 1e-9, f"tket->compactq fidelity {f}"

    # compactq -> tket and back
    rng = random.Random(5)
    for _ in range(8):
        n = 3
        ops = []
        for _ in range(12):
            if rng.random() < 0.5:
                ops.append(Gate("cx", (), tuple(rng.sample(range(n), 2))))
            else:
                nm = rng.choice(("h", "rz", "t", "s"))
                ops.append(Gate(nm, (rng.uniform(0, 6.3),) if nm == "rz"
                                else (), (rng.randrange(n),)))
        c = Circuit(n, ops)
        back = from_pytket(to_pytket(c))
        f = f2(to_qiskit(c), to_qiskit(back))
        assert f > 1 - 1e-9, f"bridge round-trip fidelity {f}"


ALL = [
    ("qiskit: PassManager output == direct compactq call",
     test_qiskit_pass_matches_direct_call),
    ("qiskit: proof verdict preserved as pass output",
     test_qiskit_pass_proof_metadata),
    ("qiskit: level-0 and approximate paths exact",
     test_qiskit_pass_level0_and_approximate_paths),
    ("pytket: pass output == direct compactq call",
     test_tket_pass_matches_direct_call),
    ("pytket: last_proof metadata + SequencePass integration",
     test_tket_pass_proof_metadata_and_sequence),
    ("bridges: pytket bridge exact in both directions",
     test_pytket_bridge_exact_both_directions),
]


def main() -> int:
    print(f"compactq transpiler-pass integration suite ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll integration tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
