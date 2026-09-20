"""Decision-diagram (QMDD-style) unitary prover — exact proof beyond the
dense ceiling.

The dense prover (compactq.equivalence) builds the full 2^n x 2^n matrix,
which caps exact proof at 8 qubits.  This module represents the same
matrix as a weighted decision diagram: equal submatrices collapse to
identical shared node structure, so structured circuits (QFT, adders,
Grover, Trotter, QAOA rings — the whole practical benchmark corpus) are
verified EXACTLY at widths where no dense matrix can exist.

Design (single canonical structure, tolerance-gated verdicts):

* A node at level L covers the top L+1 basis bits; its four successors
  cover the 2x2 block decomposition at bit L.  Level -1 is the scalar
  terminal.
* Normalization: each node's first nonzero successor weight (fixed scan
  order) is factored out to the incoming edge, making the structure
  canonical; equal submatrices share nodes through the unique table.
* Weights are hashed on a 1e-12 rounding grid so near-equal submatrices
  share (best-effort sharing).  The verdict is NEVER structural: it is
  the Hilbert-Schmidt fidelity |Tr(A+ B)|/d > 1 - tol — the same
  decision rule and tolerance as the dense prover, computed by a
  memoized recursion over the two diagrams.  Grid effects can therefore
  only affect speed, never soundness.
* A node budget bounds memory/time; exhausting it raises DDOverflow and
  the caller falls to the next prover in the cascade (decline, never
  guess) — the same contract as the rest of compactq.

Gate application needs only three primitives: 1q mixing (linear
combination of blocks), diagonal row-phasing (CZ/CP), and CX expressed
as H(target) . CZ . H(target) (CZ is symmetric).  SWAP is the standard
CX triple; multi-qubit gates are expanded first.
"""
from __future__ import annotations

from .circuit import Circuit

__all__ = ["DDOverflow", "dd_fidelity", "check_equivalent_dd",
           "dd_node_count", "DD_MAX_QUBITS"]

DD_MAX_QUBITS = 32          # hard stop; the cascade falls to randomized beyond
_ZERO_EPS = 1e-14           # edge weight below this is structurally zero
_GRID = 12                  # weight hashing grid (decimal places)


class DDOverflow(Exception):
    """Node budget exhausted — the diagram did not stay compact."""


def _q(w: complex) -> tuple:
    return (round(w.real, _GRID), round(w.imag, _GRID))


class _DD:
    """Hash-consed diagram store: unique table + per-pass memo tables +
    mark-sweep collection of dead intermediates (each gate rebuilds the
    diagram; without GC the table keeps every superseded node and
    memory explodes even when the live diagram stays tiny)."""

    def __init__(self, n: int, max_nodes: int = 400_000):
        self.n = n
        self.max_nodes = max_nodes
        self.unique: dict = {}
        self.nodes: list = []           # id -> (level, w00, c00, w01, c01,
        self.keys: list = []            #          w10, c10, w11, c11) + key
        self.free: list = []            # recycled ids
        self.allocated: set = set()     # currently occupying slots
        self.live_count = 0
        self.gc_roots: list = []        # (weight, id) pairs kept alive
        self.add_memo: dict = {}
        self.root = (0.0 + 0.0j, 0)

    def _check_budget(self):
        if len(self.allocated) > 8 * self.max_nodes:
            raise DDOverflow(
                f"decision diagram exceeded {8 * self.max_nodes} nodes")

    def sweep(self, roots):
        """Free every node unreachable from `roots` (terminal id 0 is
        implicit).  Called between gates, never mid-construction.  Stale
        add-memo entries die with the sweep: they may reference ids a
        later mk reuses for a different node."""
        mark = set()
        stack = [nid for (_, nid) in roots if nid]
        while stack:
            nid = stack.pop()
            if nid in mark:
                continue
            mark.add(nid)
            nd = self.nodes[nid - 1]
            for cid in (nd[2], nd[4], nd[6], nd[8]):
                if cid and cid not in mark:
                    stack.append(cid)
        for nid in list(self.allocated):
            if nid not in mark:
                i = nid - 1
                self.unique.pop(self.keys[i], None)
                self.nodes[i] = None
                self.free.append(nid)
        self.allocated = mark
        self.live_count = len(mark)
        self.add_memo = {}

    # -- structure ---------------------------------------------------------
    def mk(self, level: int, ch):
        """Make (or reuse) the node for four weighted children; returns
        the weighted pair (normalization_weight, id)."""
        w00, c00, w01, c01, w10, c10, w11, c11 = ch
        first = None
        for w in (w00, w01, w10, w11):
            if abs(w) > _ZERO_EPS:
                first = w
                break
        if first is None:
            return (0.0 + 0.0j, 0)
        inv = 1.0 / first
        k00 = w00 * inv
        k01 = w01 * inv
        k10 = w10 * inv
        k11 = w11 * inv
        k00 = k00 if abs(k00) > _ZERO_EPS else 0.0 + 0.0j
        k01 = k01 if abs(k01) > _ZERO_EPS else 0.0 + 0.0j
        k10 = k10 if abs(k10) > _ZERO_EPS else 0.0 + 0.0j
        k11 = k11 if abs(k11) > _ZERO_EPS else 0.0 + 0.0j
        c00 = c00 if k00 != 0 else 0
        c01 = c01 if k01 != 0 else 0
        c10 = c10 if k10 != 0 else 0
        c11 = c11 if k11 != 0 else 0
        key = (level, _q(k00), c00, _q(k01), c01, _q(k10), c10, _q(k11), c11)
        nid = self.unique.get(key)
        if nid is None:
            self._check_budget()
            if self.free:
                nid = self.free.pop()
                self.nodes[nid - 1] = (level, k00, c00, k01, c01,
                                       k10, c10, k11, c11)
                self.keys[nid - 1] = key
            else:
                nid = len(self.nodes) + 1        # 0 is the terminal
                self.nodes.append((level, k00, c00, k01, c01,
                                   k10, c10, k11, c11))
                self.keys.append(key)
            self.allocated.add(nid)
            self.unique[key] = nid
        return (first, nid)

    def get(self, nid: int):
        if nid == 0:
            return (-1, 0.0 + 0.0j, 0, 0.0 + 0.0j, 0,
                    0.0 + 0.0j, 0, 0.0 + 0.0j, 0)
        return self.nodes[nid - 1]

    # -- primitives ----------------------------------------------------------
    def add(self, a, b):
        """Weighted sum of two same-level weighted diagrams.  Memoized on
        (ids, exact weights) — hash-consing makes identical terms recur
        with identical float weights, so hits are the norm."""
        wa, ia = a
        wb, ib = b
        if abs(wa) <= _ZERO_EPS:
            return b
        if abs(wb) <= _ZERO_EPS:
            return a
        if ia == ib:                          # (wa + wb) . A
            return (wa + wb, ia)
        la = self.get(ia)[0]
        if la < 0:                            # two terminals
            return (wa + wb, 0)
        key = (ia, ib, wa, wb)
        hit = self.add_memo.get(key)
        if hit is not None:
            return hit
        ga = self.get(ia)
        gb = self.get(ib)
        ch = []
        for s in range(4):
            ch.append(self.add((wa * ga[1 + 2 * s], ga[2 + 2 * s]),
                               (wb * gb[1 + 2 * s], gb[2 + 2 * s])))
        flat = [x for pair in ch for x in pair]
        res = self.mk(la, flat)
        self.add_memo[key] = res
        return res

    def scale(self, w, a):
        if abs(w) <= _ZERO_EPS:
            return (0.0 + 0.0j, 0)
        return (a[0] * w, a[1])

    def apply_1q(self, a, wire: int, g, memo: dict):
        """Apply the 2x2 gate g=(m00, m01, m10, m11) on `wire`."""
        w, nid = a
        if abs(w) <= _ZERO_EPS:
            return (0.0 + 0.0j, 0)
        hit = memo.get(nid)
        if hit is not None:
            return (hit[0] * w, hit[1])
        nd = self.get(nid)
        level = nd[0]
        orig = ((nd[1], nd[2]), (nd[3], nd[4]), (nd[5], nd[6]),
                (nd[7], nd[8]))
        if level == wire:
            m00, m01, m10, m11 = g
            n00 = self.add(self.scale(m00, orig[0]), self.scale(m01, orig[2]))
            n01 = self.add(self.scale(m00, orig[1]), self.scale(m01, orig[3]))
            n10 = self.add(self.scale(m10, orig[0]), self.scale(m11, orig[2]))
            n11 = self.add(self.scale(m10, orig[1]), self.scale(m11, orig[3]))
        else:                                  # level > wire
            n00 = self.apply_1q(orig[0], wire, g, memo)
            n01 = self.apply_1q(orig[1], wire, g, memo)
            n10 = self.apply_1q(orig[2], wire, g, memo)
            n11 = self.apply_1q(orig[3], wire, g, memo)
        if (n00, n01, n10, n11) == orig:
            res = (1.0 + 0.0j, nid)            # unchanged subtree: keep it
        else:
            res = self.mk(level, [x for pair in (n00, n01, n10, n11)
                                  for x in pair])
        memo[nid] = res
        return (res[0] * w, res[1])

    def _neg_rows(self, a, wire: int, phase: complex, memo: dict):
        """Multiply the rows where bit `wire` is 1 by `phase`."""
        w, nid = a
        if abs(w) <= _ZERO_EPS:
            return (0.0 + 0.0j, 0)
        key = (nid, wire)
        hit = memo.get(key)
        if hit is not None:
            return (hit[0] * w, hit[1])
        nd = self.get(nid)
        level = nd[0]
        orig = ((nd[1], nd[2]), (nd[3], nd[4]), (nd[5], nd[6]),
                (nd[7], nd[8]))
        if level == wire:
            n00, n01 = orig[0], orig[1]
            n10 = self.scale(phase, orig[2])
            n11 = self.scale(phase, orig[3])
        else:
            # level > wire: rows where bit `wire` is 1 live in ALL FOUR
            # blocks (each block's rows split again by the lower bits)
            n00 = self._neg_rows(orig[0], wire, phase, memo)
            n01 = self._neg_rows(orig[1], wire, phase, memo)
            n10 = self._neg_rows(orig[2], wire, phase, memo)
            n11 = self._neg_rows(orig[3], wire, phase, memo)
        if (n00, n01, n10, n11) == orig:
            res = (1.0 + 0.0j, nid)            # unchanged subtree: keep it
        else:
            res = self.mk(level, [x for pair in (n00, n01, n10, n11)
                                  for x in pair])
        memo[key] = res
        return (res[0] * w, res[1])

    def apply_diag2(self, a, hi: int, lo: int, phase: complex):
        """Apply diag(1, 1, 1, phase) on wires {hi, lo}."""
        if hi == lo:
            raise ValueError("diagonal 2q gate needs distinct wires")
        if hi < lo:
            hi, lo = lo, hi
        w, nid = a
        if abs(w) <= _ZERO_EPS:
            return (0.0 + 0.0j, 0)
        memo: dict = {}
        rec_memo: dict = {}

        def rec(pnid: int):
            hit = rec_memo.get(pnid)
            if hit is not None:
                return hit
            nd = self.get(pnid)
            lv = nd[0]
            orig = ((nd[1], nd[2]), (nd[3], nd[4]), (nd[5], nd[6]),
                    (nd[7], nd[8]))
            if lv == hi:
                n00, n01 = orig[0], orig[1]
                n10 = self._neg_rows(orig[2], lo, phase, memo)
                n11 = self._neg_rows(orig[3], lo, phase, memo)
            else:
                # lv > hi: bit-hi=1 rows live in ALL FOUR blocks, so the
                # diag2 treatment recurses into every slot (edge weights
                # ride along; zero children stay terminal)
                def sub(pair):
                    wt, cid = pair
                    if abs(wt) <= _ZERO_EPS:
                        return pair
                    r = rec(cid)
                    return (r[0] * wt, r[1])
                n00, n01, n10, n11 = (sub(orig[0]), sub(orig[1]),
                                      sub(orig[2]), sub(orig[3]))
            if (n00, n01, n10, n11) == orig:
                res = (1.0 + 0.0j, pnid)     # unchanged subtree: keep it
            else:
                res = self.mk(lv, [x for pair in (n00, n01, n10, n11)
                                   for x in pair])
            rec_memo[pnid] = res
            return res

        res = rec(nid)
        return (res[0] * w, res[1])

    def apply_cx(self, a, c: int, t: int):
        """CX(c, t) = H(t) . CZ(c, t) . H(t)  (fresh memo per pass: the
        diagram changes between the two H applications)."""
        h = (complex(0.0, 0.7071067811865476),   # i/√2
             complex(0.0, 0.7071067811865476),   # i/√2
             complex(0.0, 0.7071067811865476),   # i/√2
             complex(0.0, -0.7071067811865476))  # -i/√2  (= i . H)
        a = self.apply_1q(a, t, h, {})
        a = self.apply_diag2(a, max(c, t), min(c, t), -1.0 + 0.0j)
        a = self.apply_1q(a, t, h, {})
        return a

    def apply_swap(self, a, x: int, y: int):
        """SWAP via the exact CX triple."""
        for c, t in ((x, y), (y, x), (x, y)):
            a = self.apply_cx(a, c, t)
        return a

    # -- equivalence ---------------------------------------------------------
    def inner(self, a, b, memo: dict):
        """|Tr(A+ B)| contribution: Tr(A+ B) for two weighted diagrams."""
        wa, ia = a
        wb, ib = b
        if abs(wa) <= _ZERO_EPS or abs(wb) <= _ZERO_EPS:
            return 0.0 + 0.0j
        key = (ia, ib)
        hit = memo.get(key)
        if hit is not None:
            return hit[0] * _cconj(wa) * wb
        la = self.get(ia)[0]
        if la < 0:                             # both terminals
            struct = 1.0 + 0.0j
        else:
            ga = self.get(ia)
            gb = self.get(ib)
            struct = 0.0 + 0.0j
            for s in range(4):
                wA = ga[1 + 2 * s]
                cA = ga[2 + 2 * s]
                wB = gb[1 + 2 * s]
                cB = gb[2 + 2 * s]
                if abs(wA) <= _ZERO_EPS or abs(wB) <= _ZERO_EPS:
                    continue
                # the recursive call already carries this slot's edge
                # weights conj(wA)*wB — they are the slot's contribution
                struct += self.inner((wA, cA), (wB, cB), memo)
        memo[key] = (struct, 0)
        return struct * _cconj(wa) * wb


def _cconj(z: complex) -> complex:
    return z.conjugate()


def _identity(store: _DD, n: int):
    ident = (1.0 + 0.0j, 0)
    for level in range(n):
        ident = store.mk(level, (1.0 + 0.0j, ident[1], 0.0 + 0.0j, 0,
                                 0.0 + 0.0j, 0, 1.0 + 0.0j, ident[1]))
    return ident


def build_dd(circ: Circuit, max_nodes: int = 400_000) -> _DD:
    """Build `circ`'s unitary as a decision diagram.  Raises ValueError
    for unsupported gates, DDOverflow past the node budget."""
    from .linalg import gate_matrix
    from .mcx import expand_gate

    n = circ.num_qubits
    if n < 1:
        raise ValueError("dd prover needs >= 1 qubit")
    if n > DD_MAX_QUBITS:
        raise ValueError(f"dd prover limited to {DD_MAX_QUBITS} qubits")
    dd = _DD(n, max_nodes)
    ident = _identity(dd, n)
    dd.gc_roots = [ident]
    dd.root = _apply_all(dd, ident, circ, gate_matrix, expand_gate)
    dd.gc_roots = [dd.root]
    dd.sweep([dd.root])
    return dd


def _apply_all(dd, cur, circ, gate_matrix, expand_gate):
    for g in circ.ops:
        cur = _apply_gate(dd, cur, g, gate_matrix, expand_gate)
        # collect superseded intermediates lazily: sweeping only when the
        # dead frontier grows past the live diagram keeps the bookkeeping
        # cost off small gates while bounding memory to ~4x live nodes
        if len(dd.allocated) > 4 * dd.live_count + 1024:
            dd.sweep(list(dd.gc_roots) + [cur])
            if dd.live_count > dd.max_nodes:
                raise DDOverflow(
                    f"decision diagram exceeded {dd.max_nodes} live nodes")
    return cur


def _apply_gate(dd, cur, g, gate_matrix, expand_gate):
    if len(g.qubits) == 1:
        return dd.apply_1q(cur, g.qubits[0], gate_matrix(g), {})
    if len(g.qubits) == 2:
        name = g.name
        if name == "cx":
            return dd.apply_cx(cur, g.qubits[0], g.qubits[1])
        if name == "cz":
            return dd.apply_diag2(cur, max(g.qubits), min(g.qubits),
                                  -1.0 + 0.0j)
        if name == "cp":
            import math
            lam = float(g.params[0]) if g.params else 0.0
            return dd.apply_diag2(cur, max(g.qubits), min(g.qubits),
                                  complex(math.cos(lam), math.sin(lam)))
        if name == "swap":
            return dd.apply_swap(cur, g.qubits[0], g.qubits[1])
        if name in ("ecr", "iswap"):
            from .native import NATIVE_EXPANSION
            for gg in NATIVE_EXPANSION[name]:
                cur = _apply_gate(dd, cur, gg, gate_matrix, expand_gate)
            return cur
        raise ValueError(f"dd prover: unsupported gate {name!r}")
    ex = expand_gate(g)
    if ex is None:
        raise ValueError(f"dd prover: unsupported gate {g.name!r}")
    for gg in ex:
        cur = _apply_gate(dd, cur, gg, gate_matrix, expand_gate)
    return cur


def dd_fidelity(circ_a: Circuit, circ_b: Circuit,
                max_nodes: int = 400_000) -> float:
    """|Tr(A+ B)|/d between the two circuits, computed on decision
    diagrams sharing one node store (common submatrices collapse, which
    makes the inner product cheap).  Raises DDOverflow / ValueError —
    the caller declines."""
    from .linalg import gate_matrix
    from .mcx import expand_gate

    if circ_a.num_qubits != circ_b.num_qubits:
        raise ValueError("qubit count mismatch")
    n = circ_a.num_qubits
    if n < 1:
        raise ValueError("dd prover needs >= 1 qubit")
    if n > DD_MAX_QUBITS:
        raise ValueError(f"dd prover limited to {DD_MAX_QUBITS} qubits")
    store = _DD(n, max_nodes)
    ident = _identity(store, n)
    store.gc_roots = [ident]
    root_a = _apply_all(store, ident, circ_a, gate_matrix, expand_gate)
    # keep A alive while B is built in the same store (shared nodes make
    # the inner product cheap; A's diagram must not be collected)
    store.gc_roots = [root_a]
    root_b = _apply_all(store, _identity(store, n), circ_b,
                        gate_matrix, expand_gate)
    store.gc_roots = [root_a, root_b]
    store.sweep(store.gc_roots)
    val = store.inner(root_a, root_b, {})
    d = 1 << n
    return abs(val) / d


def check_equivalent_dd(circ_a: Circuit, circ_b: Circuit, tol: float = 1e-7,
                        max_nodes: int = 400_000) -> bool:
    """True iff the circuits are equal up to global phase, proven on
    decision diagrams within tolerance `tol`."""
    return dd_fidelity(circ_a, circ_b, max_nodes) > 1.0 - tol


def dd_node_count(circ: Circuit, max_nodes: int = 400_000) -> int:
    """Live diagram nodes for `circ`'s unitary (tests / reporting)."""
    store = build_dd(circ, max_nodes)
    return store.live_count
