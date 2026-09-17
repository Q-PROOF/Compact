"""Deterministic search-guided optimization.

Tries several exact optimization strategies and keeps the best verified
result.  No randomness — candidate pipelines are fixed and every candidate is
re-verified against the input when the qubit count allows it.
"""
from __future__ import annotations

from .circuit import Circuit
from .equivalence import check_equivalent
from .optimize import optimize, optimize_deep, u3_fold
from .transforms import commute_cancel, peephole, slide_1q, swap_template


def _score(c: Circuit):
    return (c.two_qubit_count(), len(c.ops), c.depth())


def _deep(c: Circuit, rounds: int) -> Circuit:
    cur = c
    for _ in range(rounds):
        nxt = peephole(swap_template(slide_1q(commute_cancel(cur))))
        if (len(nxt.ops), nxt.depth()) >= (len(cur.ops), cur.depth()):
            break
        cur = nxt
    return cur




def _circuit_features(circ: Circuit) -> dict:
    """Cheap structural features for pass-selection."""
    n = circ.num_qubits
    counts = {}
    diag = 0
    cx = 0
    clifford = 0
    for g in circ.ops:
        counts[g.name] = counts.get(g.name, 0) + 1
        if len(g.qubits) == 2:
            cx += 1
        elif g.name in ("rz", "p", "s", "sdg", "z", "t", "tdg"):
            diag += 1
        if g.name in ("h", "x", "y", "z", "s", "sdg", "cx", "cz", "swap",
                      "sx", "sxdg"):
            clifford += 1
    total = max(1, len(circ.ops))
    return {
        "qubits": n,
        "gates": len(circ.ops),
        "cx_density": cx / total,
        "diagonal_density": diag / total,
        "clifford_fraction": clifford / total,
    }


def _select_candidates(feats: dict, verify: bool) -> set:
    """Decision table: which search candidates to run for this circuit kind.

    Calibrated on QASMBench (see benchmarks): the full candidate set is kept
    for small dense circuits; big/sparse/Clifford-dominated circuits skip
    the expensive low-yield candidates.
    """
    sel = {"template", "parity"}
    if verify:
        sel.add("clifford")
    if feats["gates"] <= 120:
        sel.update({"crosspair", "clifford", "parity"})
    if feats["clifford_fraction"] > 0.8:
        sel.update({"clifford"})
    if feats["cx_density"] > 0.25 and feats["gates"] <= 800:
        sel.add("crosspair")
    return sel



def optimize_search(circ: Circuit, verify: bool | None = None, depth: int = 2,
                    fidelity_tolerance: float = 1.0,
                    max_gates: int = 4000) -> Circuit:
    """Optimize with several deterministic pipelines and keep the smallest.

    Every candidate is produced by exact passes; the winner is re-verified
    against the input when `verify` allows, otherwise the plain `optimize`
    result is returned.  Guarantees: never wrong, never larger than plain
    `optimize`.

    fidelity_tolerance < 1 enables approximate mode: two-qubit blocks may be
    re-synthesized into a lower-CX class whose average gate fidelity
    (|Tr(U_ref^dag U)|/d) is at least `fidelity_tolerance`.  The exact
    whole-circuit proof is then impossible by definition, so verification is
    restricted to the exact baseline candidates and the per-block fidelity
    guarantee; the total circuit infidelity is bounded by the number of
    approximated blocks times (1 - fidelity_tolerance).  An exact baseline
    is always computed first, and an approximate candidate wins only when it
    is strictly smaller.
    """
    from .equivalence import _MAX_QUBITS as _PROOF_QUBITS
    approx = fidelity_tolerance < 1.0
    if verify is None:
        verify = circ.num_qubits <= _PROOF_QUBITS and not approx
    # scale guardrail: very large circuits get the cheap pipeline only
    if len(circ.ops) > max_gates:
        return optimize(circ, verify=False)

    base = optimize(circ, verify=verify)
    deep = optimize_deep(circ, verify=verify,
                         fid_tol=1 - (1 - fidelity_tolerance) if approx else 1 - 1e-13)
    if _score(deep) < _score(base):
        base = deep
    best, best_score = base, _score(base)

    # candidate 1: deep multi-round pass composition (no full re-verify here;
    # the whole candidate is verified below)
    deep = _deep(circ, rounds=3 * max(1, depth))
    if _score(deep) < best_score:
        if verify and not check_equivalent(circ, deep):
            return best
        best, best_score = deep, _score(deep)

    # NOTE: an older revision had a "reversed-gate-order" candidate here
    # (optimize the reversed circuit, then reverse back).  It was UNSOUND:
    # reversing a gate list does not preserve the operator, so whenever
    # optimize changed the decomposition the candidate silently became a
    # different unitary.  Removed - see tests (verify=False exactness).

    feats = _circuit_features(circ)
    sel = _select_candidates(feats, verify)

    # candidate 3: template reorderings feeding the exact pass chain
    from .templates import template_pass
    cand = template_pass(best)
    cand = optimize(cand, verify=False)
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)

    # candidate 3b: phase-polynomial re-synthesis of diagonal cores
    # (each replacement is fidelity-verified against the window unitary).
    if "parity" in sel:
        from .parity import parity_pass
        cand = parity_pass(best)
    else:
        cand = best
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)

    # candidate 4: cross-pair commutative merge + KAK on the merged blocks
    from .kak import kak_pass, merge_crosspair
    merged = merge_crosspair(best) if "crosspair" in sel else best
    cand = optimize(merged, verify=False)
    cand = kak_pass(cand, force=False)
    cand = optimize(cand, verify=False)
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)

    # candidate 4b: unitary-verified window merge + KAK
    from .winmerge import merge_by_unitary
    merged2 = merge_by_unitary(best)
    if merged2 is not best:
        cand = optimize(merged2, verify=False)
        cand = kak_pass(cand, force=False)
        cand = optimize(cand, verify=False)
        if _score(cand) < best_score:
            if not verify or check_equivalent(circ, cand):
                best, best_score = cand, _score(cand)

    # candidate 5: Clifford resynthesis (tableau-proven exact at any qubit
    # count - blocks are verified with clifford_equal inside the pass).
    if "clifford" in sel:
        from .clifford import clifford_optimize
        cand = clifford_optimize(best)
    else:
        cand = best
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)

    # candidate 5b: Clifford-structure reveal - rewrite Clifford-valued u3/rz
    # runs into canonical H/S/X, then Clifford-resynthesize (tableau-proven).
    # Helps Clifford-dense inputs; never accepted unless strictly better.
    from .clifford import cliffordize, clifford_pass as _cpass
    cand = _cpass(cliffordize(best))
    cand = optimize(cand, verify=False)
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)
    # final polish: fold 1q runs into single-u3 form (exact, gate-count only).
    # Runs LAST because u3 gates are opaque to the commutation passes.
    cand = u3_fold(best)
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)
    if _score(cand) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand)

    return best
