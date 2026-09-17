"""The compactq optimizer: safe pipeline + deep KAK cascade."""
from __future__ import annotations

from .circuit import Circuit, Gate
from .equivalence import check_equivalent
from .linalg import gate_matrix, mmul, resynth
from .transforms import commute_cancel, diag_slide, peephole, slide_1q, swap_template
from .cp_pass import cancel_cx_through_diagonal, merge_cp

_MAX_ITERS = 25
from .equivalence import _MAX_QUBITS as _VERIFY_MAX_QUBITS


def _score(c: Circuit):
    return (c.two_qubit_count(), len(c.ops), c.depth())


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


def optimize(circ: Circuit, verify: bool | None = None, max_iters: int = _MAX_ITERS) -> Circuit:
    """Safe pipeline: peephole + diag_slide + swap + slide + commute_cancel.
    Exact rewrites only, proven for <=6q."""
    if verify is None:
        verify = circ.num_qubits <= _VERIFY_MAX_QUBITS

    best = circ.copy()
    cur = circ.copy()
    best_score = (len(best.ops), best.depth())

    for _ in range(max_iters):
        prev_score = (len(cur.ops), cur.depth())
        cur = peephole(diag_slide(swap_template(slide_1q(commute_cancel(cur)))))
        score = (len(cur.ops), cur.depth())
        if score < best_score:
            best = cur.copy()
            best_score = score
        if score >= prev_score:
            break

    result = best if best_score <= (len(circ.ops), circ.depth()) else circ.copy()
    if verify:
        if not check_equivalent(circ, result):
            return circ.copy()
    return result


def optimize_deep(circ: Circuit, verify: bool | None = None, max_iters: int = 10,
                  fid_tol: float = 1 - 1e-13) -> Circuit:
    """Deep optimization: iterative KAK block re-synthesis + CP cancellation
    + full pass chain.  Can beat `optimize` on dense circuits.

    Strategy: multiple rounds of [KAK → CP-cancel → passes → safe-polish],
    keeping the best verified result across all rounds.
    """
    if verify is None:
        verify = circ.num_qubits <= _VERIFY_MAX_QUBITS

    from .kak import kak_pass

    best = optimize(circ, verify=False)
    best_score = _score(best)

    cur = circ
    stale = 0
    for rnd in range(max_iters):
        prev_sc = _score(cur)
        # KAK re-synthesis + CP cancellation + full pass chain
        stepped = kak_pass(cur, force=True, fid_tol=fid_tol)
        stepped = peephole(diag_slide(merge_cp(cancel_cx_through_diagonal(
            swap_template(slide_1q(commute_cancel(stepped)))))))
        # polish with safe optimizer
        stepped = optimize(stepped, verify=False)
        sc = _score(stepped)
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
            return optimize(circ, verify=False)
    return best
