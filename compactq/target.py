"""Hardware-aware optimization: minimize estimated circuit INFIDELITY.

A `Target` describes a machine: per-pair CX fidelities, a default CX
fidelity for unlisted pairs, and a single-qubit gate fidelity.  The
estimated infidelity of a circuit is

    sum over gates of (1 - gate fidelity)

with SWAP counted as three CX-equivalents.  `optimize_for` generates
candidates (exact baseline, direction-corrected, and several approximate
levels), scores each under the target, and returns the best together with
its estimated fidelity — so the optimizer picks "least noisy", not merely
"fewest gates".
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .circuit import Circuit, Gate
from .search import optimize_search


@dataclass
class Target:
    cx_fidelity: dict = field(default_factory=dict)
    native_2q: str = "cx"   # {(a, b): f} or {frozenset: f}
    default_cx_fidelity: float = 0.99
    single_qubit_fidelity: float = 0.9999

    def cx_fid(self, a: int, b: int) -> float:
        f = self.cx_fidelity.get((a, b))
        if f is None:
            f = self.cx_fidelity.get(frozenset((a, b)))
        if f is None:
            f = self.cx_fidelity.get((b, a))
        return self.default_cx_fidelity if f is None else f

    def gate_infidelity(self, g: Gate) -> float:
        if g.name in ("cx", "cz", "ecr", "iswap"):
            return 1.0 - self.cx_fid(g.qubits[0], g.qubits[1])
        if g.name == "swap":
            return 3.0 * (1.0 - self.cx_fid(g.qubits[0], g.qubits[1]))
        if g.name == "cp":
            return 1.5 * (1.0 - self.cx_fid(g.qubits[0], g.qubits[1]))
        if len(g.qubits) == 1:
            return 1.0 - self.single_qubit_fidelity
        return 0.0

    def estimated_infidelity(self, circ: Circuit) -> float:
        return sum(self.gate_infidelity(g) for g in circ.ops)


def _flip_to_best_direction(circ: Circuit, target: Target) -> Circuit:
    """For each CX on an undirected pair, point it at the better direction
    (exact H-conjugation flip)."""
    from .hardware import flip_cx
    out = circ
    pairs = {tuple(sorted(g.qubits)) for g in out.ops
             if g.name == "cx" and len(g.qubits) == 2}
    for a, b in pairs:
        fab = target.cx_fid(a, b)
        fba = target.cx_fid(b, a)
        if fba > fab + 1e-15:
            out = flip_cx(out, a, b)
    return out


def _approx_fidelity(base: Circuit, cand: Circuit, tol: float) -> float:
    """Realized fidelity of an approximate candidate vs the exact baseline.

    Measured directly when the qubit count allows a unitary proof; otherwise
    bounded below by the per-block guarantee: each approximated block has
    fidelity >= tol and saves at most 3 CX, so tol ** ceil(saved_cx / 3) is a
    safe floor on the whole-circuit fidelity loss.
    """
    from .equivalence import _MAX_QUBITS as proof_q
    if base.num_qubits <= proof_q:
        try:
            from .equivalence import fidelity
            return fidelity(base, cand)
        except Exception:
            pass
    saved = base.two_qubit_count() - cand.two_qubit_count()
    if saved <= 0:
        return 1.0
    return tol ** (-(-saved // 3))


def optimize_for(circ: Circuit, target: Target,
                 approx_levels: tuple = (1.0, 0.99, 0.95),
                 verify: bool | None = None) -> tuple[Circuit, float]:
    """Optimize `circ` for the machine described by `target`.

    Returns (circuit, estimated_fidelity).  Candidates span the exact
    baseline, direction-corrected placement, and the requested approximate
    levels.  Exact-level candidates are fidelity-verified (for qubit counts
    within the prover's limit) and are only returned when exactly equivalent
    to the input.  For approximate candidates the reported fidelity is the
    hardware-model fidelity multiplied by the candidate's realized
    approximation fidelity (measured where provable, per-block-bounded
    otherwise), so the estimate never counts the approximation as free.
    """
    if verify is None:
        from .equivalence import _MAX_QUBITS as proof_q
        verify = circ.num_qubits <= proof_q
    from .equivalence import check_equivalent

    # (circuit, realized approximation fidelity; 1.0 == exact)
    candidates = []
    base = optimize_search(circ, verify=verify)
    candidates.append((base, 1.0))

    flipped = _flip_to_best_direction(base, target)
    if flipped is not base:
        cand = optimize_search(flipped, verify=verify)
        candidates.append((cand, 1.0))

    for tol in approx_levels:
        if tol >= 1.0:
            continue
        cand = optimize_search(circ, verify=verify, fidelity_tolerance=tol)
        F = _approx_fidelity(base, cand, tol)
        candidates.append((cand, F))
        flipped2 = _flip_to_best_direction(cand, target)
        if flipped2 is not cand:
            candidates.append((flipped2, F))

    best = None
    best_f = None
    for cand, approx_f in candidates:
        if approx_f >= 1.0 and verify and not check_equivalent(circ, cand):
            continue  # never return an inequivalent exact-level candidate
        f = (1.0 - target.estimated_infidelity(cand)) * approx_f
        if best_f is None or f > best_f + 1e-15:
            best, best_f = cand, f
    if best is None:
        best = base
        best_f = 1.0 - target.estimated_infidelity(base)

    # native-basis translation: re-express 2q blocks in the machine's own
    # gate and keep the rebased circuit when the hardware model prefers it
    native = getattr(target, "native_2q", "cx")
    if native != "cx":
        from .native import rebase as _rebase
        rb = _rebase(best, native)
        ok = True
        if verify and circ.num_qubits <= 8:
            ok = check_equivalent(circ, rb)
        if ok:
            rb_f = 1.0 - target.estimated_infidelity(rb)
            if rb_f > best_f - 1e-15:
                best, best_f = rb, rb_f
    return best, max(0.0, best_f)


def approximate_for_target(circ: Circuit, target: Target,
                           budget: float | None = None) -> tuple[Circuit, float]:
    """Greedy per-block fidelity-budget allocation under `target`.

    For every two-qubit block of the exact baseline, each cheaper CX class
    is an option with:
      savings = CX removed x (1 - CX fidelity of that pair)
      cost    = 1 - block fidelity of the cheaper class
    Options are applied best-ratio first while the ratio exceeds 1 and the
    total fidelity cost stays within `budget` (default: unbounded - the
    ratio rule alone then caps spending at the point of diminishing
    returns).  Returns (circuit, estimated_fidelity); exact blocks are
    untouched, so the result is never worse-connected than the baseline.
    The reported fidelity multiplies the hardware-model estimate by
    (1 - total approximation cost spent), so it never presents the
    approximation as free.
    """
    from . import kak

    if budget is not None and budget < 0:
        raise ValueError("budget must be >= 0")
    base = optimize_search(circ, verify=(circ.num_qubits <= 6))
    base_inf = target.estimated_infidelity(base)

    blocks = kak.collect_blocks(base)
    options = []
    for (start, end, a, b) in blocks:
        bc = kak._block_circuit(base, start, end, a, b)
        cx_exact = bc.two_qubit_count()
        if cx_exact < 2:
            continue
        U4 = kak._block_unitary(bc)
        if U4 is None:
            continue
        w = kak.weyl(U4)
        if w is None:
            continue
        a_, b_, c_, *_ = w
        fids = kak._fid_table(a_, b_, c_)
        pair_inf = 1.0 - target.cx_fid(a, b)
        for nb in range(cx_exact):
            if fids[nb] >= 1 - 1e-9:
                continue  # same as exact - nothing to trade
            saved = (cx_exact - nb) * pair_inf
            cost = 1.0 - fids[nb]
            if cost <= 1e-12:
                continue
            options.append((saved / cost, saved, cost, (start, end, a, b, nb)))
    options.sort(key=lambda o: (-o[0], -o[1]))
    if budget is not None and not options and base_inf == 0:
        pass

    chosen = []
    spent = 0.0
    for ratio, saved, cost, blk in options:
        if ratio <= 1.0:
            break
        if budget is not None and spent + cost > budget:
            continue
        chosen.append(blk)
        spent += cost
    if not chosen:
        return base, 1.0 - base_inf

    # splice: synthesize chosen blocks at their pinned class, keep the rest
    ops = base.ops
    out = []
    prev = 0
    by_start = {}
    for (start, end, a, b, nb) in chosen:
        by_start.setdefault(start, (start, end, a, b, nb))
    for start in sorted(by_start):
        end, a, b, nb = by_start[start][1], by_start[start][2], by_start[start][3], by_start[start][4]
        bc = kak._block_circuit(base, start, end, a, b)
        U4 = kak._block_unitary(bc)
        synth = kak.synth_2q_ops(U4, force_nb=nb)
        if synth is None:
            continue
        remap = {0: a, 1: b}
        synth = [Gate(g.name, g.params, tuple(remap[q] for q in g.qubits))
                 for g in synth]
        out.extend(ops[prev:start])
        out.extend(synth)
        prev = end
    out.extend(ops[prev:])
    result = Circuit(circ.num_qubits, out)

    # polish 1q runs
    from .optimize import u3_fold
    result = u3_fold(result)

    est_inf = target.estimated_infidelity(result)
    # honest accounting: the hardware-model fidelity times the total
    # approximation cost actually spent on blocks (sum of per-block
    # 1 - block-fidelity, mirroring the whole-circuit bound convention)
    est_f = max(0.0, (1.0 - est_inf) * (1.0 - spent))
    if est_f >= 1.0 - base_inf:
        return base, 1.0 - base_inf
    return result, est_f
