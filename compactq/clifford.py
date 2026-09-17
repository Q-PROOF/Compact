"""Clifford block resynthesis (Aaronson-Gottesman synthesis).

Ports the AG tableau-reduction synthesis (Aaronson & Gottesman 2004; the
same construction as Qiskit's `synth_clifford_ag`): a working copy of the
block's tableau is reduced to the identity by applying H/S/CX/SWAP/X/Z
gates, and the emitted circuit is the inverse of that reduction sequence.

Correctness contract:
- the reduction is self-checked (the working tableau must be in the exact
  identity pattern up to row phases when reduction completes, else the pass
  refuses to fire);
- every synthesized candidate is verified with `clifford_equal` — an exact
  proof at ANY qubit count (no 2^n cost), so Clifford optimization is
  honestly verified well beyond the 6-qubit unitary limit.
"""
from __future__ import annotations

from .circuit import Circuit, Gate
from .stabilizer import (
    Tableau,
    _expand_cz_swap,
    _gate_1q_table,
    clifford_equal,
    tableau,
)

_INV = {"h": "h", "s": "sdg", "sdg": "s", "x": "x", "z": "z",
        "cx": "cx", "swap": "swap", "sx": "sxdg", "sxdg": "sx"}


def _is_clifford_1q(g: Gate, cache={}):
    key = (g.name, g.params)
    if key not in cache:
        cache[key] = _gate_1q_table(g)
    return cache[key]


class _Reducer:
    """Applies gates to a working tableau while recording them."""

    def __init__(self, tab: Tableau):
        self.t = tab
        self.red: list[Gate] = []

    def cx(self, c, t):
        self.t.apply_cx(c, t)
        self.red.append(Gate("cx", (), (c, t)))

    def swap(self, a, b):
        self.t.apply_swap(a, b)
        self.red.append(Gate("swap", (), (a, b)))

    def _1q(self, name, q):
        table = _is_clifford_1q(Gate(name, (), (q,)))
        self.t.apply_1q_table(q, table)
        self.red.append(Gate(name, (), (q,)))

    def h(self, q):
        self._1q("h", q)

    def s(self, q):
        self._1q("s", q)

    def z(self, q):
        self._1q("z", q)

    def x(self, q):
        self._1q("x", q)


def _synth_gates_from_tableau(tab: Tableau) -> list[Gate] | None:
    """AG reduction: gates g with tableau(circuit(g)) == tab (up to phase)."""
    n = tab.n
    t = tab.copy()
    r = _Reducer(t)

    def destab_x(q):
        return t.x[q]

    def destab_z(q):
        return t.z[q]

    def stab_x(q):
        return t.x[n + q]

    def stab_z(q):
        return t.z[n + q]

    def set_qubit_x_true(q):
        if t.x[q][q]:
            return
        for i in range(q + 1, n):
            if destab_x(q)[i]:
                r.swap(i, q)
                return
        # no non-zero X element: apply a Hadamard to bring a Z over
        for i in range(q, n):
            if destab_z(q)[i]:
                r.h(i)
                if i != q:
                    r.swap(i, q)
                return

    def set_row_x_zero(q):
        row_x = destab_x(q)
        row_z = destab_z(q)
        for i in range(q + 1, n):
            if row_x[i]:
                r.cx(q, i)
        if any(row_z[j] for j in range(q, n)):
            if not row_z[q]:
                r.s(q)
            for i in range(q + 1, n):
                if row_z[i]:
                    r.cx(i, q)
            r.s(q)

    def set_row_z_zero(q):
        row_x = stab_x(q)
        row_z = stab_z(q)
        if any(row_z[j] for j in range(q + 1, n)):
            for i in range(q + 1, n):
                if row_z[i]:
                    r.cx(i, q)
        if any(row_x[j] for j in range(q, n)):
            r.h(q)
            for i in range(q + 1, n):
                if row_x[i]:
                    r.cx(q, i)
            if row_z[q]:
                r.s(q)
            r.h(q)

    for q in range(n):
        set_qubit_x_true(q)
        set_row_x_zero(q)
        set_row_z_zero(q)

    # leftover row signs: after the pattern checks each row is a single-qubit
    # X_q (destab) or Z_q (stab) form; Z(q) / X(q) anticommute with exactly
    # those, flipping their -1 signs back to +1.
    for q in range(n):
        if t.e[q] % 4:
            r.z(q)
        if t.e[n + q] % 4:
            r.x(q)

    # self-check: the working tableau must now match the identity pattern
    # exactly (the GF(2) corrector below is a second safety net)
    for q in range(n):
        for j in range(n):
            if t.x[q][j] != (1 if j == q else 0):
                return None
            if t.z[q][j] != 0:
                return None
            if t.x[n + q][j] != 0:
                return None
            if t.z[n + q][j] != (1 if j == q else 0):
                return None
        if t.e[q] % 4 or t.e[n + q] % 4:
            return None

    # the synthesis is the inverse of the reduction sequence
    out = [Gate(_INV[g.name], g.params, g.qubits) for g in reversed(r.red)]

    # sign corrector: whatever residual -1 signs the emitted circuit carries
    # (AG phase bookkeeping is convention-sensitive), solve for appended
    # Z(q)/X(q) gates that flip exactly the differing rows.  The symplectic
    # tableau is full rank, so a solution always exists; if anything is
    # off-pattern we return None and the caller keeps the original block.
    try:
        t_out = tableau(Circuit(n, out))
    except Exception:
        return None
    corr = _sign_corrector(tab, t_out)
    if corr is None:
        return None
    out.extend(corr)
    return out


def _sign_corrector(tab_in: Tableau, tab_out: Tableau):
    """Gates Z(q)/X(q) whose appended conjugation flips exactly the row signs
    where tab_out and tab_in differ by -1.

    Appending Z(q) flips every row with an X factor on q (column q of the
    x-matrix); appending X(q) flips rows with a Z factor on q.  The tableau
    is symplectic full-rank, so [X | Z] is invertible over GF(2) and a
    solution always exists.  Returns None when the difference is not
    sign-only (pattern mismatch).
    """
    n = tab_in.n
    d = []
    for i in range(2 * n):
        diff = (tab_out.e[i] - tab_in.e[i]) % 4
        if diff % 2:
            return None  # odd phase difference: not a sign-level mismatch
        d.append(diff // 2)  # 1 iff the row sign differs (-1)
    if not any(d):
        return []
    for i in range(2 * n):
        if tab_out.x[i] != tab_in.x[i] or tab_out.z[i] != tab_in.z[i]:
            return None
    rows = 2 * n
    M = [[0] * (2 * n) for _ in range(rows)]
    for i in range(rows):
        for q in range(n):
            M[i][q] = tab_in.x[i][q]
            M[i][n + q] = tab_in.z[i][q]
    aug = [M[i][:] + [d[i]] for i in range(rows)]
    pivots = []
    row = 0
    for col in range(2 * n):
        piv = next((r2 for r2 in range(row, rows) if aug[r2][col]), None)
        if piv is None:
            continue
        aug[row], aug[piv] = aug[piv], aug[row]
        for r2 in range(rows):
            if r2 != row and aug[r2][col]:
                aug[r2] = [(a ^ b) for a, b in zip(aug[r2], aug[row])]
        pivots.append((row, col))
        row += 1
    for r2 in range(row, rows):
        if aug[r2][-1]:
            return None  # inconsistent (cannot happen for symplectic rows)
    u = [0] * (2 * n)
    for (r2, col) in pivots:
        u[col] = aug[r2][-1]
    gates = []
    for q in range(n):
        if u[q]:
            gates.append(Gate("z", (), (q,)))
    for q in range(n):
        if u[n + q]:
            gates.append(Gate("x", (), (q,)))
    return gates


# ---------------------------------------------------------------------------
# block collection + pass
# ---------------------------------------------------------------------------

def collect_clifford_blocks(circ: Circuit, min_cx: int = 3):
    """Maximal contiguous windows that are entirely Clifford.

    Returns (start, end, tableau) tuples; the tableau is the block's
    simualtion and is reused for synthesis.  cz expands internally; swap and
    all Clifford 1q gates are native.  Blocks with fewer than `min_cx`
    CX-equivalents are not reported (nothing to win on 2q count).
    """
    from .circuit import Gate as _G

    ops = []
    for g in circ.ops:  # expand cz -> h cx h once, up front
        if g.name == "cz":
            a, b = g.qubits
            ops.extend([_G("h", (), (b,)), _G("cx", (), (a, b)), _G("h", (), (b,))])
        else:
            ops.append(g)

    blocks = []
    i = 0
    n_ops = len(ops)
    while i < n_ops:
        g = ops[i]
        ok_start = (
            (len(g.qubits) == 1 and _is_clifford_1q(g) is not None)
            or g.name == "cx"
            or g.name == "swap"
        )
        if not ok_start:
            i += 1
            continue
        tab = Tableau(circ.num_qubits)
        j = i
        n_cx = 0
        while j < n_ops:
            g = ops[j]
            if len(g.qubits) == 1:
                table = _is_clifford_1q(g)
                if table is None:
                    break
                tab.apply_1q_table(g.qubits[0], table)
            elif g.name == "cx":
                tab.apply_cx(g.qubits[0], g.qubits[1])
                n_cx += 1
            elif g.name == "swap":
                tab.apply_swap(g.qubits[0], g.qubits[1])
                n_cx += 1  # swap costs like a CX in metrics
            else:
                break
            j += 1
        if n_cx >= min_cx and j - i >= 3:
            blocks.append((i, j, tab, n_cx))
            i = j
        else:
            i = j if j > i else i + 1
    return blocks


def clifford_pass(circ: Circuit, force: bool = False) -> Circuit:
    """Resynthesize maximal Clifford blocks into canonical H/S/CX form.

    Every candidate is proven exact with `clifford_equal` (a tableau proof,
    valid at any qubit count).  Blocks are only replaced when the result is
    not worse on the (2q, gates) score, unless force=True.
    """
    blocks = collect_clifford_blocks(circ)
    if not blocks:
        return circ
    ops = circ.ops
    out = []
    prev = 0
    changed = False
    proven = True  # every accepted replacement carries a tableau proof
    for (start, end, tab, n_cx) in blocks:
        out.extend(ops[prev:start])
        prev = end
        synth = _synth_gates_from_tableau(tab)
        if synth is None:
            out.extend(ops[start:end])
            continue
        block_circ = Circuit(circ.num_qubits, ops[start:end])
        cand = Circuit(circ.num_qubits, synth)
        try:
            if not clifford_equal(block_circ, cand):
                out.extend(ops[start:end])
                continue
        except ValueError:
            out.extend(ops[start:end])
            continue
        block_2q = block_circ.two_qubit_count()
        cand_2q = cand.two_qubit_count()
        if force:
            accept = True
        else:
            accept = (cand_2q < block_2q
                      or (cand_2q == block_2q and len(synth) < end - start))
        if accept:
            out.extend(synth)
            changed = True
        else:
            out.extend(ops[start:end])
    out.extend(ops[prev:])
    if not changed:
        return circ
    return Circuit(circ.num_qubits, out)


def clifford_optimize(circ: Circuit) -> Circuit:
    """Clifford resynthesis + safe-pipeline polish (exact)."""
    from .transforms import commute_cancel, peephole, slide_1q, swap_template

    stepped = clifford_pass(circ, force=True)
    stepped = peephole(swap_template(slide_1q(commute_cancel(stepped))))
    # a second short Clifford round can fire after the polish merged runs
    stepped = clifford_pass(stepped, force=True)
    stepped = peephole(swap_template(slide_1q(commute_cancel(stepped))))
    return stepped


_CLIFFORD_GENS = None


def cliffordize(circ: Circuit) -> Circuit:
    """Re-express Clifford-valued 1q gates (u3/rz/p at Clifford angles) as
    canonical H/S/X products so block collection can see the runs.
    Non-Clifford gates pass through untouched.  Every rewrite is proven
    by matrix equality up to global phase."""
    from .linalg import gate_matrix, same_up_to_phase
    global _CLIFFORD_GENS
    if _CLIFFORD_GENS is None:
        gens = []
        for nm in ("h", "s", "x"):
            gens.append((nm, (), gate_matrix(Gate(nm, (), (0,)))))
        seen = {}
        frontier = [((), ())]
        # BFS over products of generators (mod global phase), depth <= 4
        elems = {(): ((), ())}
        prods = [((), ())]
        found = {}
        for depth in range(5):
            nxt = []
            for seq, mat in prods:
                for nm, prm, gm in gens:
                    mm = _mm2(gm, mat) if mat else gm
                    key = _phase_key2(mm)
                    if key in found:
                        continue
                    fseq = (seq + (nm,))  # circuit order: left-to-right
                    found[key] = fseq
                    nxt.append((fseq, mm))
            prods = nxt
        _CLIFFORD_GENS = found
    out = []
    changed = 0
    for g in circ.ops:
        if len(g.qubits) != 1 or g.name in ("h", "s", "sdg", "x", "y", "z",
                                              "sdg", "id"):
            out.append(g)
            continue
        m = None
        try:
            m = gate_matrix(g)
        except Exception:
            m = None
        if m is None:
            out.append(g)
            continue
        key = _phase_key2(m)
        seq = _CLIFFORD_GENS.get(key)
        if seq is None:
            out.append(g)
            continue
        for nm in seq:
            out.append(Gate(nm, (), g.qubits))
        changed += 1
    return Circuit(circ.num_qubits, out)


def _phase_key(m):
    """Phase-insensitive canonical key of a 2x2 unitary (or None)."""
    import cmath
    if not m or m[0] == 0:
        return None
    ph = m[0] / abs(m[0])
    return tuple(round(x, 9) for x in
                 (v / ph for v in m)) if all(abs(v) < 1e9 for v in m) else None

def _phase_key2(m):
    """Phase-insensitive canonical key of a 2x2 unitary matrix."""
    import cmath
    m00 = complex(m[0])
    if abs(m00) < 1e-12:
        m00 = complex(m[1])
        if abs(m00) < 1e-12:
            return None
    ph = m00 / abs(m00)
    return tuple(complex(round((complex(v) / ph).real, 9),
                         round((complex(v) / ph).imag, 9)) for v in m)


def _mm2(a, b):
    """2x2 matrix product for flat 4-tuples (a,b,c,d)."""
    a0, a1, a2, a3 = complex(a[0]), complex(a[1]), complex(a[2]), complex(a[3])
    b0, b1, b2, b3 = complex(b[0]), complex(b[1]), complex(b[2]), complex(b[3])
    return (a0 * b0 + a1 * b2, a0 * b1 + a1 * b3,
            a2 * b0 + a3 * b2, a2 * b1 + a3 * b3)
