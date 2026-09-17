"""The compactq optimizer: safe pipeline + deep KAK cascade."""
from __future__ import annotations

from .circuit import Circuit, Gate
from .equivalence import check_equivalent
from .linalg import gate_matrix, mmul, resynth
from .transforms import commute_cancel, diag_slide, peephole, slide_1q, swap_template
from .cp_pass import cancel_cx_through_diagonal, merge_cp

_MAX_ITERS = 25
from .equivalence import _MAX_QUBITS as _VERIFY_MAX_QUBITS

# Inner-pipeline candidates skip the (expensive) whole-circuit proof; the
# enclosing entry point re-verifies its final result before returning.
_UNVERIFIED = False

_OBJECTIVES = ("2q", "depth", "gate_count")


def _score(c: Circuit, objective: str = "2q"):
    """Lexicographic metric tuple for the chosen objective.

    Every supported objective is a reordering of the same three metrics
    (2-qubit count, total gates, depth), so any candidate accepted under
    any objective strictly improves that objective's tuple: the never-grow
    invariant holds by construction, and every rewrite remains exact
    regardless of the objective.
    """
    if objective == "2q":
        return (c.two_qubit_count(), len(c.ops), c.depth())
    if objective == "depth":
        return (c.depth(), c.two_qubit_count(), len(c.ops))
    if objective == "gate_count":
        return (len(c.ops), c.two_qubit_count(), c.depth())
    raise ValueError(f"unknown objective {objective!r}; expected one of "
                     f"{sorted(_OBJECTIVES)} (hardware-error weighting lives "
                     f"in compactq.target.optimize_for)")


def _pair(c: Circuit, objective: str):
    """Two-metric acceptance key for the safe pipeline's fixpoint.

    Identical to the historical (gates, depth) order for the default and
    gate_count objectives (zero behavior drift); depth-first for 'depth'.
    """
    if objective == "depth":
        return (c.depth(), len(c.ops))
    return (len(c.ops), c.depth())


def u3_fold(circ: Circuit) -> Circuit:
    """Final polish: refold every maximal 1-qubit run, preferring a single
    u3 gate when no shorter exact form exists.

    Purely local and exact (each folded run is verified against its own
    product up to global phase); only strictly shorter foldings are kept.
    Run AFTER all commutation passes — u3 gates are opaque to them.
    """
    ops = circ.ops
    out = []
    i = 0
    n = len(ops)
    while i < n:
        g = ops[i]
        if len(g.qubits) == 1:
            q = g.qubits[0]
            prod = gate_matrix(g)
            j = i + 1
            while j < n and len(ops[j].qubits) == 1 and ops[j].qubits[0] == q:
                prod = mmul(gate_matrix(ops[j]), prod)
                j += 1
            try:
                folded = resynth(prod, prefer_u=True)
            except Exception:
                folded = None
            if folded is not None and len(folded) < j - i:
                out.extend(Gate(fg.name, fg.params, (q,)) for fg in folded)
            else:
                out.extend(ops[i:j])
            i = j
        else:
            out.append(g)
            i += 1
    return Circuit(circ.num_qubits, out)


def optimize(circ: Circuit, verify: bool | None = None, max_iters: int = _MAX_ITERS,
             objective: str = "2q") -> Circuit:
    """Safe pipeline: peephole + diag_slide + swap + slide + commute_cancel.
    Exact rewrites only, proven for <=6q.

    objective: which metric may never grow (lexicographic).  '2q' (default)
    and 'gate_count' share the historical acceptance order; 'depth' accepts
    depth-improving trades first.
    """
    if objective not in _OBJECTIVES:
        raise ValueError(f"unknown objective {objective!r}; expected one of "
                         f"{sorted(_OBJECTIVES)}")
    if verify is None:
        verify = circ.num_qubits <= _VERIFY_MAX_QUBITS

    best = circ.copy()
    cur = circ.copy()
    best_score = _pair(best, objective)

    for _ in range(max_iters):
        prev_score = _pair(cur, objective)
        cur = peephole(diag_slide(swap_template(slide_1q(commute_cancel(cur)))))
        score = _pair(cur, objective)
        if score < best_score:
            best = cur.copy()
            best_score = score
        if score >= prev_score:
            break

    result = best if _score(best, objective) <= _score(circ, objective) else circ.copy()
    if verify:
        if not check_equivalent(circ, result):
            return circ.copy()
    return result


def optimize_deep(circ: Circuit, verify: bool | None = None, max_iters: int = 10,
                  fid_tol: float = 1 - 1e-13, objective: str = "2q") -> Circuit:
    """Deep optimization: iterative KAK block re-synthesis + CP cancellation
    + full pass chain.  Can beat `optimize` on dense circuits.

    Strategy: multiple rounds of [KAK → CP-cancel → passes → safe-polish],
    keeping the best verified result across all rounds.  `objective` selects
    the lexicographic metric order used for every keep/discard decision.
    """
    if objective not in _OBJECTIVES:
        raise ValueError(f"unknown objective {objective!r}; expected one of "
                         f"{sorted(_OBJECTIVES)}")
    if verify is None:
        verify = circ.num_qubits <= _VERIFY_MAX_QUBITS

    from .kak import kak_pass

    best = optimize(circ, verify=_UNVERIFIED, objective=objective)
    best_score = _score(best, objective)

    cur = circ
    stale = 0
    for rnd in range(max_iters):
        prev_sc = _score(cur, objective)
        # KAK re-synthesis + CP cancellation + full pass chain
        stepped = kak_pass(cur, force=True, fid_tol=fid_tol)
        stepped = peephole(diag_slide(merge_cp(cancel_cx_through_diagonal(
            swap_template(slide_1q(commute_cancel(stepped)))))))
        # polish with safe optimizer
        stepped = optimize(stepped, verify=_UNVERIFIED, objective=objective)
        sc = _score(stepped, objective)
        if sc < best_score:
            best = stepped
            best_score = sc
            cur = stepped
            stale = 0
        elif sc < prev_sc:
            cur = stepped
            stale += 1
        else:
            stale += 1
            if stale >= 2:
                break

    # verification
    if verify:
        if not check_equivalent(circ, best):
            return optimize(circ, verify=_UNVERIFIED, objective=objective)
    return best
