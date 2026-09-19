"""Verification-API tests: compactq.verify() tiers and evidence dict.

Runs standalone (`python tests/test_verify_api.py`, zero dependencies) and
is wired into CI alongside the other suites.
"""
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


def test_verify_equivalent_pairs():
    import compactq
    from compactq.benchmarks import qft, ghz
    # Clifford pair -> tier 3, algebraic
    r = compactq.verify(ghz(5), ghz(5))
    assert r["equivalent"] is True and r["tier"] == 3, r
    assert r["method"] == "clifford_tableau" and r["global_phase_ignored"] is True
    # general pair -> tier 2, dense unitary
    r = compactq.verify(qft(4), compactq.optimize(qft(4)))
    assert r["equivalent"] is True and r["tier"] == 2, r
    assert r["method"] == "full_unitary"
    assert isinstance(r["runtime_ms"], int)


def test_verify_inequivalent_pair():
    import compactq
    from compactq.benchmarks import qft, ghz
    r = compactq.verify(qft(4), ghz(4))
    assert r["equivalent"] is False and r["tier"] == 2, r


def test_verify_rejections_and_shape():
    import compactq
    from compactq.benchmarks import qft
    r = compactq.verify(qft(4), qft(3))
    assert r["equivalent"] is False and r["method"] == "shape_mismatch"
    # unknown-objective rejection is part of the objective API contract
    try:
        compactq.optimize(qft(3), objective="bogus")
        raise AssertionError("expected ValueError")
    except ValueError:
        pass


def test_verify_tier_coverage_documented():
    import compactq
    tiers = compactq.VERIFICATION_TIERS
    assert set(tiers) == {0, 1, 2, 3, 4, 5}
    # guarantee vocabulary: exact-proven / verified-randomized / unverified
    assert tiers[0] == "unverified"
    assert tiers[1] == "verified-randomized"
    assert tiers[2].startswith("exact-proven")
    assert tiers[3].startswith("exact-proven")
    # T4 is SHIPPED (compositional) — it must not be marked roadmap
    assert tiers[4] == "exact-proven (compositional certificates)"
    assert "roadmap" in tiers[5]


def test_proof_status_tier_mapping():
    import compactq
    from compactq.verify import proof_status_tier
    assert proof_status_tier("exact-unitary") == 2
    assert proof_status_tier("compositional") == 4
    assert proof_status_tier("randomized-exact") == 1
    assert proof_status_tier("unverified") == 0
    assert proof_status_tier("approximate") == 0
    assert proof_status_tier("unknown-status") == 0


ALL = [
    ("verify(): equivalent pairs report correct tiers", test_verify_equivalent_pairs),
    ("verify(): inequivalent pair detected (tier 2)", test_verify_inequivalent_pair),
    ("verify(): shape mismatch + objective rejection", test_verify_rejections_and_shape),
    ("verify(): tier table T0-T5 documented", test_verify_tier_coverage_documented),
]


def main() -> int:
    print(f"compactq verification-API tests ({len(ALL)} tests)")
    for name, fn in ALL:
        check(name, fn)
    if FAILED:
        print(f"\n{len(FAILED)} FAILED: {FAILED}")
        return 1
    print("\nAll verification-API tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
