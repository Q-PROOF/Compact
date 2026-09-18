"""Trust-layer tests: T4 compositional verification, certificate format,
verify CLI subcommand, and coupling-map presets.

Runs standalone (`python tests/test_trust_layer.py`, zero dependencies)
and is wired into CI alongside the other suites.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

FAILED = []


def check(name, fn):
    try:
        fn()
        print(f"  PASS  {name}")
    except Exception as e:
        FAILED.append(name)
        print(f"  FAIL  {name}: {type(e).__name__}: {e}")


def _blocks_100q():
    from compactq import Circuit, Gate
    a = [Gate("h", (), (0,))] + [Gate("cx", (), (j, j + 1)) for j in range(3)]
    b = [Gate("ry", (0.7,), (50,)), Gate("cx", (), (50, 51)),
         Gate("rz", (0.3,), (52,)), Gate("cx", (), (51, 52))]
    return Circuit(100, a + b), a, b


def test_t4_equivalent_100q():
    import compactq
    orig, a, b = _blocks_100q()
    r = compactq.verify(orig, orig)
    assert r["equivalent"] is True and r["tier"] == 4, r
    assert r["method"] == "compositional_disjoint_blocks" and r["blocks"] == 2


def test_t4_tampered_block_named():
    import compactq
    from compactq import Circuit, Gate
    orig, a, b = _blocks_100q()
    bad = Circuit(100, a + [Gate("cx", (), (50, 51)),
                            Gate("rz", (0.3,), (52,)),
                            Gate("cx", (), (51, 52))])
    r = compactq.verify(orig, bad)
    assert r["equivalent"] is False and r["tier"] == 4, r
    assert r["failing_block_qubits"] == [50, 51, 52], r


def test_certificate_structure():
    import compactq
    from compactq.benchmarks import qft
    orig = qft(4)
    opt = compactq.optimize(orig)
    cert = compactq.build_certificate(orig, opt)
    for key in ("compiler", "version", "input_hash", "output_hash",
                "input_qubits", "optimization", "verification"):
        assert key in cert, key
    assert cert["compiler"] == "compactq"
    assert cert["verification"]["equivalent"] is True
    assert cert["optimization"]["gates_after"] <= cert["optimization"]["gates_before"]
    assert cert["input_hash"] != cert["output_hash"] or orig is opt
    # deterministic: same circuits -> same hashes
    cert2 = compactq.build_certificate(orig, opt)
    assert cert["input_hash"] == cert2["input_hash"]


def test_coupling_presets():
    import compactq
    line = compactq.coupling_preset("line", 6)
    assert len(line) == 5 and {0, 1} in line
    a2a = compactq.coupling_preset("all_to_all", 5)
    assert len(a2a) == 10
    hh = compactq.coupling_preset("heavy_hex", 16)
    assert len(hh) >= 12 and all(len(p) == 2 for p in hh)
    grid = compactq.coupling_preset("grid", 9)
    assert len(grid) == 12  # 3x3 grid: 6 horizontal + 6 vertical
    try:
        compactq.coupling_preset("nonsense", 8)
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_verify_cli_subcommand(tmp: Path):
    import compactq
    from compactq.benchmarks import qft
    from compactq.io_qasm import to_qasm
    from compactq.__main__ import main
    f1 = tmp / "orig.qasm"
    f2 = tmp / "opt.qasm"
    f3 = tmp / "cert.json"
    orig = qft(4)
    opt = compactq.optimize(orig)
    f1.write_text(to_qasm(orig), encoding="utf-8")
    f2.write_text(to_qasm(opt), encoding="utf-8")
    rc = main(["verify", str(f1), str(f2), "-o", str(f3)])
    assert rc == 0
    cert = json.loads(f3.read_text(encoding="utf-8"))
    assert cert["verification"]["equivalent"] is True
    assert cert["compiler"] == "compactq"


ALL = [
    ("T4 compositional: 100Q disjoint blocks equivalent (tier 4)",
     test_t4_equivalent_100q),
    ("T4 compositional: tampered block detected and named",
     test_t4_tampered_block_named),
    ("certificate: structure, hashes, metrics", test_certificate_structure),
    ("coupling presets: line/all-to-all/heavy-hex/grid",
     test_coupling_presets),
    ("verify CLI subcommand end-to-end", lambda: test_verify_cli_subcommand(
        Path(__file__).resolve().parent / "_trust_tmp")),
]


def main() -> int:
    tmp = Path(__file__).resolve().parent / "_trust_tmp"
    tmp.mkdir(exist_ok=True)
    print(f"compactq trust-layer tests ({len(ALL)} tests)")
    for name, fn in ALL:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as e:
            FAILED.append(name)
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll trust-layer tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
