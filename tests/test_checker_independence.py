"""Checker-independence tests (T1.1.5).

Property: random circuits → optimize_with_certificate → compactq-check
says VALID.  Mutations: every corruption of the certificate or the
optimized circuit is caught.  The checker module never imports compactq.

Standalone (`python tests/test_checker_independence.py`).
"""
import copy
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]
                       / "compactq-check" / "src"))

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


# ------------------------------------------------------------------ tests
def test_checker_never_imports_compactq():
    import sys
    for m in list(sys.modules):
        if m == "compactq" or m.startswith("compactq."):
            del sys.modules[m]
    import compactq_check.checker  # noqa: F401
    bad = [m for m in sys.modules
           if m == "compactq" or m.startswith("compactq.")]
    assert not bad, f"checker imported compactq: {bad}"


def test_property_certificates_valid():
    import compactq
    import compactq_check.checker as checker
    from compactq.cert import optimize_with_certificate
    from compactq.benchmarks import qft, ghz, random_circuit

    poly_ops = []
    for j in range(5):
        poly_ops.append(compactq.Gate("cx", (), (j, j + 1)))
        poly_ops.append(compactq.Gate("rz", (0.7 + 0.1 * j,), (j + 1,)))
    poly = compactq.Circuit(6, poly_ops)

    circuits = [qft(4), qft(6), ghz(5), random_circuit(4, 24, seed=11),
                random_circuit(6, 30, seed=12), ghz(8), poly,
                compactq.optimize(poly)]
    kinds = {}
    for circ in circuits:
        res = optimize_with_certificate(circ)
        assert res.certificate is not None, \
            f"no cert for {circ.num_qubits}q: {res.reason}"
        cert = res.certificate
        kind = cert["witness"]["kind"]
        kinds[kind] = kinds.get(kind, 0) + 1
        in_q = compactq.to_qasm(circ)
        out_q = compactq.to_qasm(res.circuit)
        v = checker.check_certificate(cert, in_q, out_q)
        assert v["verdict"] == "VALID", (circ.num_qubits, v)
    print(f"    witness kinds used: {kinds}")


def test_mutations_all_killed():
    import compactq
    import compactq_check.checker as checker
    from compactq.cert import optimize_with_certificate
    from compactq.benchmarks import ghz, qft

    base = ghz(3)
    res = optimize_with_certificate(base)
    cert = res.certificate
    assert cert is not None
    in_q = compactq.to_qasm(base)
    out_q = compactq.to_qasm(res.circuit)

    def flip_hash(d, key):
        d = copy.deepcopy(d)
        h = d[key]["sha256"]
        d[key]["sha256"] = h[:-1] + ("0" if h[-1] != "0" else "1")
        return d

    def with_n(d, n):
        d = copy.deepcopy(d)
        d["n_qubits"] = n
        return d

    def with_version(d, v):
        d = copy.deepcopy(d)
        d["cert_version"] = v
        return d

    def with_kind(d, k):
        d = copy.deepcopy(d)
        d["witness"]["kind"] = k
        return d

    def with_composition(d, comp):
        d = copy.deepcopy(d)
        d["composition"] = comp
        return d

    def without_witness(d):
        d = copy.deepcopy(d)
        d.pop("witness")
        return d

    def with_input_hash_flipped(d):
        d = copy.deepcopy(d)
        h = d["input"]["sha256"]
        d["input"]["sha256"] = ("0" if h[0] != "0" else "1") + h[1:]
        return d

    def with_output_hash_flipped(d):
        d = copy.deepcopy(d)
        h = d["output"]["sha256"]
        d["output"]["sha256"] = ("0" if h[0] != "0" else "1") + h[1:]
        return d

    mutants = [
        ("input hash altered", with_input_hash_flipped(cert), in_q, out_q),
        ("output hash altered", with_output_hash_flipped(cert), in_q, out_q),
        ("n_qubits inflated", with_n(cert, cert["n_qubits"] + 1), in_q, out_q),
        ("cert_version mutated", with_version(cert, 99), in_q, out_q),
        ("witness kind bogus", with_kind(cert, "bogus_kind"), in_q, out_q),
        ("composition mutated", with_composition(cert, "tensor_local_v1"),
         in_q, out_q),
        ("witness removed", without_witness(cert), in_q, out_q),
    ]

    # circuit-level mutants: a wrong circuit must never validate
    def drop_last_gate(t):
        lines = t.splitlines(keepends=True)
        for i in range(len(lines) - 1, -1, -1):
            s = lines[i].strip()
            if s and not s.startswith(("OPENQASM", "include", "qreg",
                                       "creg")):
                return "".join(lines[:i] + lines[i + 1:])
        return t

    def flip_last_cx(t):
        lines = t.splitlines()
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].strip().startswith("cx "):
                a, b = lines[i].strip().rstrip(";").split()[1:]
                lines[i] = f"cx {b}, {a};"
                return "\n".join(lines) + "\n"
        return t

    out_drop = drop_last_gate(out_q)
    mutants.append(("last gate dropped", cert, in_q, out_drop))
    out_flip = flip_last_cx(out_q)
    if out_flip is not None:
        mutants.append(("cx direction flipped", cert, in_q, out_flip))

    survived = []
    for label, c, i, o in mutants:
        v = checker.check_certificate(c, i, o)
        if v["verdict"] == "VALID":
            survived.append(label)
    assert not survived, f"mutants survived: {survived}"
    print(f"    mutation-kill rate: 100% ({len(mutants)} mutants)")


def test_informational_claims_ignored():
    import compactq
    import compactq_check.checker as checker
    from compactq.cert import optimize_with_certificate
    from compactq.benchmarks import ghz
    import copy

    base = ghz(3)
    res = optimize_with_certificate(base)
    cert = copy.deepcopy(res.certificate)
    cert["claims"]["exact"] = False  # tampering with claims is ignored:
    # claims are informational; the recomputed evidence is what matters
    in_q = compactq.to_qasm(base)
    out_q = compactq.to_qasm(res.circuit)
    v = checker.check_certificate(cert, in_q, out_q)
    assert v["verdict"] == "VALID", v


ALL = [
    ("checker never imports compactq", test_checker_never_imports_compactq),
    ("property: certificates checker-VALID across witness kinds",
     test_property_certificates_valid),
    ("mutation kill: every corruption caught (100%)",
     test_mutations_all_killed),
    ("informational claims field ignored safely",
     test_informational_claims_ignored),
]


def main() -> int:
    print(f"compactq checker-independence tests ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll checker-independence tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
