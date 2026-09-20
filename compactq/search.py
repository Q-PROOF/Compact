"""Deterministic search-guided optimization.

Tries several exact optimization strategies and keeps the best verified
result.  No randomness — candidate pipelines are fixed and every candidate is
re-verified against the input when the qubit count allows it.

objectives: '2q' (default), 'depth', 'gate_count' reorder the lexicographic
metric tuple used for every keep/discard decision (the chosen primary metric
can never grow); 'latency' is a depth alias; 'weighted' ranks candidates by
1.0*2q + 0.1*depth + 0.02*gates (the weighted cost never exceeds the input's).
"""
from __future__ import annotations

from .circuit import Circuit
from .equivalence import check_equivalent
from .optimize import optimize, optimize_deep, u3_fold, _UNVERIFIED, _OBJECTIVES
from .transforms import commute_cancel, peephole, slide_1q, swap_template

_SEARCH_OBJECTIVES = _OBJECTIVES + ("weighted",)

# Heavy-circuit guardrail (measured on the scalability gauntlet): beyond
# this gate count the expensive synthesis candidates (phase-polynomial
# re-synthesis, cross-pair KAK, unitary window merge, Cliffordize) can
# dominate runtime for marginal gain on wide structured circuits, so they
# are skipped and the cheap exact pipeline + template/fold candidates run.
# Small benchmark circuits (QASMBench small max ~500 gates) are unaffected.
_HEAVY_GATE_LIMIT = 1200

# Small-circuit fast path: below this gate count the candidate portfolio
# rarely pays for itself (measured on the latency bench) — the safe
# pipeline + deep round + u3 fold are run, the heavy synthesis
# candidates are skipped to keep interactive latency low.
_FAST_PATH_GATES = 20

# weighted-mode cost: NISQ-oriented — the 2-qubit count dominates, depth and
# total gates are secondary penalties.  Hardware-error weighting (per-pair
# fidelities) lives in compactq.target.optimize_for.
_WEIGHTS = {"2q": 1.0, "depth": 0.1, "gates": 0.02}


def _score(c: Circuit, objective: str = "2q"):
    if objective == "weighted":
        return (_WEIGHTS["2q"] * c.two_qubit_count()
                + _WEIGHTS["depth"] * c.depth()
                + _WEIGHTS["gates"] * len(c.ops),)
    from .optimize import _score as _opt_score
    return _opt_score(c, objective)


def _normalize(objective: str) -> str:
    if objective == "latency":  # depth is the latency proxy
        return "depth"
    if objective not in _SEARCH_OBJECTIVES:
        raise ValueError(
            f"unknown objective {objective!r}; expected one of "
            f"{sorted(_SEARCH_OBJECTIVES)} (hardware-error weighting lives "
            f"in compactq.target.optimize_for)")
    return objective


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
                    max_gates: int = 4000, objective: str = "2q") -> Circuit:
    """Optimize with several deterministic pipelines and keep the smallest.

    Every candidate is produced by exact passes; the winner is re-verified
    against the input when `verify` allows, otherwise the plain `optimize`
    result is returned.  Guarantees: never wrong, never larger than plain
    `optimize`.

    objective selects the acceptance order: '2q' (default), 'depth',
    'gate_count' are lexicographic reorderings of (2q, gates, depth) — the
    chosen primary metric never grows; 'latency' aliases 'depth'; 'weighted'
    ranks candidates by 1.0*2q + 0.1*depth + 0.02*gates and guarantees the
    weighted cost never exceeds the input's.

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
    objective = _normalize(objective)
    pipe_obj = objective if objective in _OBJECTIVES else "2q"
    from .equivalence import _MAX_QUBITS as _PROOF_QUBITS
    approx = fidelity_tolerance < 1.0
    if verify is None:
        verify = circ.num_qubits <= _PROOF_QUBITS and not approx
    # scale guardrail: very large circuits get the cheap pipeline only
    if len(circ.ops) > max_gates:
        return optimize(circ, verify=_UNVERIFIED, objective=pipe_obj)

    base = optimize(circ, verify=verify, objective=pipe_obj)
    deep = optimize_deep(circ, verify=verify, objective=pipe_obj,
                         fid_tol=1 - (1 - fidelity_tolerance) if approx else 1 - 1e-13)
    if _score(deep, objective) < _score(base, objective):
        base = deep
    best, best_score = base, _score(base, objective)

    # candidate 1: deep multi-round pass composition (no full re-verify here;
    # the whole candidate is verified below)
    deep = _deep(circ, rounds=3 * max(1, depth))
    if _score(deep, objective) < best_score:
        if verify and not check_equivalent(circ, deep):
            return best
        best, best_score = deep, _score(deep, objective)

    # small-circuit fast path: interactive latency beats marginal
    # synthesis candidates — one u3 fold, then done
    if len(circ.ops) <= _FAST_PATH_GATES:
        cand = u3_fold(best)
        if _score(cand, objective) < best_score:
            if not verify or check_equivalent(circ, cand):
                best, best_score = cand, _score(cand, objective)
        return best

    # NOTE: an older revision had a "reversed-gate-order" candidate here
    # (optimize the reversed circuit, then reverse back).  It was UNSOUND:
    # reversing a gate list does not preserve the operator, so whenever
    # optimize changed the decomposition the candidate silently became a
    # different unitary.  Removed - see tests (no-verify exactness).

    feats = _circuit_features(circ)
    sel = _select_candidates(feats, verify)
    heavy = len(circ.ops) > _HEAVY_GATE_LIMIT
    if heavy:
        sel -= {"parity", "crosspair", "clifford"}

    # candidate 3: template reorderings feeding the exact pass chain
    from .templates import template_pass
    cand = template_pass(best)
    cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    # candidate 3b: phase-polynomial re-synthesis of diagonal cores
    # (each replacement is fidelity-verified against the window unitary).
    if "parity" in sel and not heavy:
        from .parity import parity_pass
        cand = parity_pass(best)
    else:
        cand = best
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    # candidate 4: cross-pair commutative merge + KAK on the merged blocks
    from .kak import kak_pass, merge_crosspair
    merged = merge_crosspair(best) if "crosspair" in sel else best
    cand = optimize(merged, verify=_UNVERIFIED, objective=pipe_obj)
    cand = kak_pass(cand, force=False)
    cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    # candidate 4b: unitary-verified window merge + KAK
    from .winmerge import merge_by_unitary
    merged2 = merge_by_unitary(best) if not heavy else best
    if merged2 is not best:
        cand = optimize(merged2, verify=_UNVERIFIED, objective=pipe_obj)
        cand = kak_pass(cand, force=False)
        cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
        if _score(cand, objective) < best_score:
            if not verify or check_equivalent(circ, cand):
                best, best_score = cand, _score(cand, objective)

    # candidate 5: Clifford resynthesis (tableau-proven exact at any qubit
    # count - blocks are verified with clifford_equal inside the pass).
    if "clifford" in sel:
        from .clifford import clifford_optimize
        cand = clifford_optimize(best)
    else:
        cand = best
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    # candidate 5b: Clifford-structure reveal - rewrite Clifford-valued u3/rz
    # runs into canonical H/S/X, then Clifford-resynthesize (tableau-proven).
    # Helps Clifford-dense inputs; never accepted unless strictly better.
    from .clifford import cliffordize, clifford_pass as _cpass
    cand = _cpass(cliffordize(best))
    cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    # candidate 6: permutation-aware KAK - per-block free SWAP orientation
    # with the residual permutation paid once as SWAPs (SWAP elision, the
    # Qiskit-L3 trick).  Verified per block inside the pass.
    if not heavy:
        from .permkak import permutation_kak_pass
        pk = permutation_kak_pass(best)
        if pk is not None:
            cand = optimize(pk, verify=_UNVERIFIED, objective=pipe_obj)
            if _score(cand, objective) < best_score:
                if not verify or check_equivalent(circ, cand):
                    best, best_score = cand, _score(cand, objective)

    # candidate 6b: pair-locality packing - exact reordering (disjoint-wire
    # gates commute) that packs same-pair gates into contiguous windows so
    # KAK resynthesis sees bigger blocks.  Verified whole-circuit below.
    if not heavy:
        from .pairpack import pair_pack
        packed = pair_pack(best)
        if len(packed.ops) and [ (g.name, g.qubits) for g in packed.ops ] != \
               [ (g.name, g.qubits) for g in best.ops ]:
            cand = optimize(packed, verify=_UNVERIFIED, objective=pipe_obj)
            cand = kak_pass(cand, force=False)
            cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
            if _score(cand, objective) < best_score:
                if not verify or check_equivalent(circ, cand):
                    best, best_score = cand, _score(cand, objective)

    # candidate 7: Clifford+T normal form - factors the circuit into
    # Clifford layers (tableau-exact) + shared parity phase layers;
    # whole-circuit-verified inside the pass and only accepted when the
    # 2q count strictly drops.
    if not heavy and feats.get("clifford_fraction", 0.0) > 0.3:
        from .cliffordt import cliffordt_pass
        cand = cliffordt_pass(best)
        if cand is not None:
            cand = optimize(cand, verify=_UNVERIFIED, objective=pipe_obj)
            if _score(cand, objective) < best_score:
                if not verify or check_equivalent(circ, cand):
                    best, best_score = cand, _score(cand, objective)

    # final polish: fold 1q runs into single-u3 form (exact, gate-count only).
    # Runs LAST because u3 gates are opaque to the commutation passes.
    cand = u3_fold(best)
    if _score(cand, objective) < best_score:
        if not verify or check_equivalent(circ, cand):
            best, best_score = cand, _score(cand, objective)

    return best
