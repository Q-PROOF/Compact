"""Gate semantics for the checker: dense matrices, Clifford canonical
forms, and phase-polynomial extraction — implemented fresh from the
standard definitions (AG 2004 tableau rules; Amy/Maslov/Mosca phase
polynomials).  No code shared with the optimizer.
"""
from __future__ import annotations

import math

I2 = [[1, 0], [0, 1]]
X = [[0, 1], [1, 0]]
Y = [[0, -1j], [1j, 0]]
Z = [[1, 0], [0, -1]]
H = [[1 / math.sqrt(2), 1 / math.sqrt(2)],
     [1 / math.sqrt(2), -1 / math.sqrt(2)]]
S = [[1, 0], [0, 1j]]
SDG = [[1, 0], [0, -1j]]
T = [[1, 0], [0, complex(2 ** -0.5, 2 ** -0.5)]]
TDG = [[1, 0], [0, complex(2 ** -0.5, -(2 ** -0.5))]]
SX = [[complex(0.5, 0.5), complex(0.5, -0.5)],
      [complex(0.5, -0.5), complex(0.5, 0.5)]]
SXDG = [[complex(0.5, -0.5), complex(0.5, 0.5)],
        [complex(0.5, 0.5), complex(0.5, -0.5)]]


def rz(t):
    return [[1, 0], [0, complex(math.cos(t), math.sin(t))]]


def rx(t):
    c, s = math.cos(t / 2), math.sin(t / 2)
    return [[c, -1j * s], [-1j * s, c]]


def ry(t):
    c, s = math.cos(t / 2), math.sin(t / 2)
    return [[c, -s], [s, c]]


def p(t):
    return [[1, 0], [0, complex(math.cos(t), math.sin(t))]]


def u3(t, phi, lam):
    c, s = math.cos(t / 2), math.sin(t / 2)
    return [[c, -s * complex(math.cos(lam), math.sin(lam))],
            [s * complex(math.cos(phi), math.sin(phi)),
             c * complex(math.cos(phi + lam), math.sin(phi + lam))]]


def u1(lam):
    return [[1, 0], [0, complex(math.cos(lam), math.sin(lam))]]


CX = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
CZ = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]]
SWAP = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
CCX = [[1.0 if (i == j or (i in (6, 7) and j in (6, 7))) else 0.0
        for j in range(8)] for i in range(8)]
CCX[6][7] = 1.0
CCX[7][6] = 1.0
CCX[6][6] = 0.0
CCX[7][7] = 0.0


def _kron(a, b):
    ra, ca = len(a), len(a[0])
    rb, cb = len(b), len(b[0])
    out = [[0j] * (ca * cb) for _ in range(ra * rb)]
    for i in range(ra):
        for j in range(ca):
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k][j * cb + l] = a[i][j] * b[k][l]
    return out


def gate_matrix(name, params):
    t = params[0] if len(params) > 0 else 0.0
    table = {
        "id": I2, "x": X, "y": Y, "z": Z, "h": H, "s": S, "sdg": SDG,
        "t": T, "tdg": TDG, "sx": SX, "sxdg": SXDG,
        "rx": rx(t), "ry": ry(t), "rz": rz(t), "p": p(t), "u1": u1(t),
    }
    if name in table:
        return table[name]
    if name in ("u3", "u"):
        t2 = params[0] if len(params) > 0 else 0.0
        phi = params[1] if len(params) > 1 else 0.0
        lam = params[2] if len(params) > 2 else 0.0
        return u3(t2, phi, lam)
    if name == "cx":
        return CX
    if name == "cz":
        return CZ
    if name == "swap":
        return SWAP
    if name == "ccx":
        return CCX
    if name in ("cp", "cu1"):
        lam = params[0] if params else 0.0
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                [0, 0, 0, complex(math.cos(lam), math.sin(lam))]]
    raise ValueError(f"unsupported gate {name!r}")


def apply_dense(U, g, ws, n):
    """Return G·U where G acts on `ws`.  Wire conventions: the GLOBAL
    basis is little-endian (wire w <-> bit w, matching compactq's QASM
    semantics), while the gate matrix's LOCAL basis is big-endian
    (ws[0] is the most significant local wire — control before target
    for CX), matching the standard gate-matrix definitions.

    new[r][c] = sum_kpat g[rpat][kpat] * U[s][c], where rpat = r's ws
    bits and s = r with its ws bits replaced by kpat.  Every (r, c) is
    recomputed — cells the gate does not mix are re-derived too, never
    carried over stale."""
    d = 1 << n
    shifts = [w for w in ws]
    L = len(ws)
    new = [[0j] * d for _ in range(d)]
    for r in range(d):
        rpat_bits = [(r >> shifts[b]) & 1 for b in range(L)]
        rpat = 0
        for b in range(L):
            rpat |= rpat_bits[b] << (L - 1 - b)
        for c in range(d):
            acc = 0j
            for kpat in range(m_local(L)):
                s = r
                for b in range(L):
                    kb = (kpat >> (L - 1 - b)) & 1
                    s = (s & ~(1 << shifts[b])) | (kb << shifts[b])
                acc += g[rpat][kpat] * U[s][c]
            new[r][c] = acc
    return new


def m_local(L):
    return 1 << L


def build_unitary(ops, n):
    d = 1 << n
    U = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
    for (name, params, qubits) in ops:
        g = gate_matrix(name, params)
        U = apply_dense(U, g, qubits, n)
    return U


def active_wires(ops):
    """Sorted wires a gate list actually touches (idle wires carry
    identity and are factored out of dense re-derivations)."""
    ws = set()
    for (_name, _params, qubits) in ops:
        ws.update(qubits)
    return sorted(ws)


def active_components(ops):
    """Split `ops` into per-component gate lists over disjoint wire sets
    (union-find over gate qubits; op order preserved within a
    component; components sorted by their sorted wire lists).  Mirrors
    compactq.compositional.active_components so the checker re-derives
    the producer's block structure independently."""
    if not ops:
        return []
    parent = {q: q for g in ops for q in g[2]}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in ops:
        r0 = find(g[2][0])
        for q in g[2][1:]:
            parent[find(q)] = r0
    groups = {}
    for g in ops:
        groups.setdefault(find(g[2][0]), []).append(g)
    comps = []
    for _root, gops in groups.items():
        ws = sorted({q for g in gops for q in g[2]})
        remap = {w: i for i, w in enumerate(ws)}
        comps.append((ws, [(g[0], g[1], tuple(remap[q] for q in g[2]))
                           for g in gops]))
    comps.sort(key=lambda t: t[0])
    return comps


def build_unitary_active(ops, n):
    """(active_wires, unitary-on-active-width) for `ops` over `n` wires:
    idle wires are dropped and the touched wires remapped compactly."""
    ws = active_wires(ops)
    if not ws:
        return [], [[1.0]]
    remap = {w: i for i, w in enumerate(ws)}
    compact = [(name, params, tuple(remap[w] for w in qubits))
               for (name, params, qubits) in ops]
    return ws, build_unitary(compact, len(ws))


def fidelity(A, B):
    d = len(A)
    acc = 0j
    for i in range(d):
        for j in range(d):
            acc += (A[i][j].conjugate()) * B[i][j]
    return abs(acc) / d


# -------------------- Clifford canonical form --------------------
# Pauli generator images as (xmask, zmask, sign) over n qubits.


CLIFFORD_1Q = {"h", "x", "y", "z", "s", "sdg"}
CLIFFORD_2Q = {"cx", "cz", "swap"}


def is_clifford_ops(ops):
    return all(o[0] in CLIFFORD_1Q or o[0] in CLIFFORD_2Q for o in ops)


def clifford_canonical(ops, n):
    """Canonical image of the 2n Pauli generators under the Clifford
    circuit.  Returns tuple of (xmask, zmask, sign) per generator, with
    sign ∈ {0,1} (Pauli = i^sign · X^x Z^z)."""
    gens = []
    for j in range(n):
        gens.append([1 << j, 0, 0])       # X_j
    for j in range(n):
        gens.append([0, 1 << j, 0])       # Z_j
    for (name, _params, qubits) in ops:
        if name == "cx":
            c, t = qubits
            for g in gens:
                xc, zc = (g[0] >> c) & 1, (g[1] >> c) & 1
                xt, zt = (g[0] >> t) & 1, (g[1] >> t) & 1
                g[0] ^= (xt << c)
                g[1] ^= (zc << t)
            continue
        if name == "cz":
            c, t = qubits
            for g in gens:
                zt = (g[1] >> t) & 1
                g[1] ^= (zt << c)
            continue
        if name == "swap":
            a, b = qubits
            for g in gens:
                xa, xb = (g[0] >> a) & 1, (g[0] >> b) & 1
                za, zb = (g[1] >> a) & 1, (g[1] >> b) & 1
                g[0] = (g[0] & ~(1 << a) & ~(1 << b)) | (xb << a) | (xa << b)
                g[1] = (g[1] & ~(1 << a) & ~(1 << b)) | (zb << a) | (za << b)
            continue
        for g in gens:
            a = (g[0] >> qubits[0]) & 1
            b = (g[1] >> qubits[0]) & 1
            if name == "h":
                g[0], g[1] = b, a
                g[2] ^= (a & b)
            elif name == "s":
                g[1] ^= a
                g[2] ^= a
            elif name == "sdg":
                g[1] ^= a
            elif name == "z":
                g[2] ^= a
            elif name == "x":
                g[2] ^= b
            elif name == "y":
                g[2] ^= (a ^ b)
    return tuple((g[0], g[1], g[2]) for g in gens)


def clifford_equal(ops_a, ops_b, n):
    """Clifford circuits equal up to global phase iff their generator
    images match in (x,z) and their signs agree up to one global bit."""
    ca = clifford_canonical(ops_a, n)
    cb = clifford_canonical(ops_b, n)
    forms_a = [(x, z) for x, z, _ in ca]
    forms_b = [(x, z) for x, z, _ in cb]
    if forms_a != forms_b:
        return False
    signs_a = [s for _, _, s in ca]
    signs_b = [s for _, _, s in cb]
    ref = signs_a[0] ^ signs_b[0]
    return all((sa ^ sb) == ref for sa, sb in zip(signs_a, signs_b))


# -------------------- phase polynomials --------------------
PHASEPOLY_OPS = {"cx", "rz", "p", "s", "sdg", "t", "tdg", "z", "id"}


def is_phasepoly_ops(ops):
    return all(o[0] in PHASEPOLY_OPS for o in ops)


def phasepoly_canonical(ops, n, tol=1e-7):
    """Canonical parity→angle table for a CNOT+diagonal circuit at any
    width.  Returns (sorted [(parity_int, angle)], ok) with angles mod 2π."""
    parities = [1 << q for q in range(n)]
    angles = {0: 0.0}
    for (name, params, qubits) in ops:
        theta = params[0] if params else 0.0
        if name == "cx":
            c, t = qubits
            parities[t] ^= parities[c]
        else:
            key = parities[qubits[0]]
            angles[key] = (angles.get(key, 0.0) + _phase_angle(name, theta))
    table = sorted((k, _mod(v)) for k, v in angles.items())
    return table, True


def _phase_angle(name, theta):
    return {"rz": theta, "p": theta, "t": math.pi / 4,
            "tdg": -math.pi / 4, "s": math.pi / 2,
            "sdg": -math.pi / 2, "z": math.pi}.get(name, 0.0)


def _mod(x):
    y = math.fmod(x, 2 * math.pi)
    if y > math.pi:
        y -= 2 * math.pi
    if y < -math.pi:
        y += 2 * math.pi
    return y


def phasepoly_equal(ops_a, ops_b, n, tol=1e-7):
    ta, _ = phasepoly_canonical(ops_a, n, tol)
    tb, _ = phasepoly_canonical(ops_b, n, tol)
    if len(ta) != len(tb):
        return False
    for (ka, va), (kb, vb) in zip(ta, tb):
        if ka != kb:
            return False
        if abs(_mod(va - vb)) > tol:
            return False
    return True


# ------------------------------------------------ decision-diagram witness
# Independent QMDD-style re-derivation for "dd" certificates: weighted
# hash-consed diagrams, first-nonzero-slot normalization, verdict from
# the memoized overlap |Tr(A+ B)|/d at the certificate tolerance.
# Written against the certificate spec; shares no producer code.

_DD_EPS = 1e-14


class _DDB:
    __slots__ = ("nodes", "table", "addmemo")

    def __init__(self):
        self.nodes = []   # id-1 -> (level, w0,c0, w1,c1, w2,c2, w3,c3)
        self.table = {}
        self.addmemo = {}

    def mk(self, level, slots):
        pivot = None
        for w in slots[0::2]:
            if abs(w) > _DD_EPS:
                pivot = w
                break
        if pivot is None:
            return (0.0 + 0.0j, 0)
        f = 1.0 / pivot
        k = [0.0 + 0.0j if abs(w * f) <= _DD_EPS else w * f
             for w in slots[0::2]]
        cs = [c if abs(w) > _DD_EPS else 0
              for w, c in zip(k, slots[1::2])]
        key = (level, round(k[0].real, 12), round(k[0].imag, 12), cs[0],
               round(k[1].real, 12), round(k[1].imag, 12), cs[1],
               round(k[2].real, 12), round(k[2].imag, 12), cs[2],
               round(k[3].real, 12), round(k[3].imag, 12), cs[3])
        nid = self.table.get(key)
        if nid is None:
            nid = len(self.nodes) + 1
            self.nodes.append((level, k[0], cs[0], k[1], cs[1],
                               k[2], cs[2], k[3], cs[3]))
            self.table[key] = nid
        return (pivot, nid)

    def add(self, a, b):
        wa, ia = a
        wb, ib = b
        if abs(wa) <= _DD_EPS:
            return b
        if abs(wb) <= _DD_EPS:
            return a
        if ia == ib:
            return (wa + wb, ia)
        na = self.nodes[ia - 1]
        nb = self.nodes[ib - 1]
        la, lb = na[0], nb[0]
        if la < 0 and lb < 0:
            return (wa + wb, 0)
        if la < 0 or lb < 0:
            raise ValueError("dd level mismatch in checker add")
        ch = []
        for s in range(4):
            ch.append(self.add((wa * na[1 + 2 * s], na[2 + 2 * s]),
                               (wb * nb[1 + 2 * s], nb[2 + 2 * s])))
        res = self.mk(la, [x for p in ch for x in p])
        return res

    def apply_1q(self, a, wire, m):
        w, nid = a
        if nid == 0 or abs(w) <= _DD_EPS:
            return (0.0 + 0.0j, 0)
        nd = self.nodes[nid - 1]
        s = [(nd[1], nd[2]), (nd[3], nd[4]), (nd[5], nd[6]), (nd[7], nd[8])]
        if nd[0] == wire:
            n00 = self.add((m[0] * s[0][0], s[0][1]),
                           (m[1] * s[2][0], s[2][1]))
            n01 = self.add((m[0] * s[1][0], s[1][1]),
                           (m[1] * s[3][0], s[3][1]))
            n10 = self.add((m[2] * s[0][0], s[0][1]),
                           (m[3] * s[2][0], s[2][1]))
            n11 = self.add((m[2] * s[1][0], s[1][1]),
                           (m[3] * s[3][0], s[3][1]))
        else:
            def down(pair):
                r = self.apply_1q((1.0 + 0.0j, pair[1]), wire, m)
                return (r[0] * pair[0], r[1])
            n00, n01 = down(s[0]), down(s[1])
            n10, n11 = down(s[2]), down(s[3])
        res = self.mk(nd[0], [x for p in (n00, n01, n10, n11) for x in p])
        return (res[0] * w, res[1])

    def apply_cz(self, a, hi, lo, phase):
        w, nid = a
        if nid == 0 or abs(w) <= _DD_EPS:
            return (0.0 + 0.0j, 0)
        nd = self.nodes[nid - 1]

        def negrows(pair):
            wt, cid = pair
            if cid == 0 or abs(wt) <= _DD_EPS:
                return pair
            m = self.nodes[cid - 1]
            if m[0] == lo:
                ch = [(m[1], m[2]), (m[3], m[4]),
                      (phase * m[5], m[6]), (phase * m[7], m[8])]
                r = self.mk(m[0], [x for p in ch for x in p])
            else:
                ch = [negrows((m[1], m[2])), negrows((m[3], m[4])),
                      negrows((m[5], m[6])), negrows((m[7], m[8]))]
                r = self.mk(m[0], [x for p in ch for x in p])
            return (r[0] * wt, r[1])

        if nd[0] == hi:
            ch = [(nd[1], nd[2]), (nd[3], nd[4]),
                  negrows((nd[5], nd[6])), negrows((nd[7], nd[8]))]
        else:
            def down(pair):
                wt, cid = pair
                if cid == 0 or abs(wt) <= _DD_EPS:
                    return pair
                r = self.apply_cz((1.0 + 0.0j, cid), hi, lo, phase)
                return (r[0] * wt, r[1])
            ch = [down((nd[1], nd[2])), down((nd[3], nd[4])),
                  down((nd[5], nd[6])), down((nd[7], nd[8]))]
        res = self.mk(nd[0], [x for p in ch for x in p])
        return (res[0] * w, res[1])

    def overlap(self, a, b, memo):
        wa, ia = a
        wb, ib = b
        if abs(wa) <= _DD_EPS or abs(wb) <= _DD_EPS:
            return 0.0 + 0.0j
        key = (ia, ib)
        hit = memo.get(key)
        if hit is not None:
            return hit * wa.conjugate() * wb
        la = self.nodes[ia - 1][0] if ia else -1
        if la < 0:
            struct = 1.0 + 0.0j
        else:
            na = self.nodes[ia - 1]
            nb = self.nodes[ib - 1]
            struct = 0.0 + 0.0j
            for s in range(4):
                wa_s, ia_s = na[1 + 2 * s], na[2 + 2 * s]
                wb_s, ib_s = nb[1 + 2 * s], nb[2 + 2 * s]
                if abs(wa_s) <= _DD_EPS or abs(wb_s) <= _DD_EPS:
                    continue
                struct += self.overlap((wa_s, ia_s), (wb_s, ib_s), memo)
        memo[key] = struct
        return struct * wa.conjugate() * wb


def _dd_ops_unitary(b, ops, n):
    cur = _ddb_identity(b, n)
    for (name, params, qubits) in ops:
        if len(qubits) == 1:
            g = gate_matrix(name, params)
            # the checker stores matrices nested 2x2; the diagram wants
            # a flat (m00, m01, m10, m11)
            if isinstance(g[0], (list, tuple)):
                m = (g[0][0], g[0][1], g[1][0], g[1][1])
            else:
                m = (g[0], g[1], g[2], g[3])
            cur = b.apply_1q(cur, qubits[0], m)
        elif len(qubits) == 2:
            if name == "cx":
                cur = _ddb_cx(b, cur, qubits[0], qubits[1])
            elif name == "cz":
                hi, lo = max(qubits), min(qubits)
                cur = b.apply_cz(cur, hi, lo, -1.0 + 0.0j)
            elif name in ("cp", "cu1"):
                hi, lo = max(qubits), min(qubits)
                lam = params[0] if params else 0.0
                cur = b.apply_cz(cur, hi, lo,
                                 complex(math.cos(lam), math.sin(lam)))
            elif name == "swap":
                for c, t in ((qubits[0], qubits[1]),
                             (qubits[1], qubits[0]),
                             (qubits[0], qubits[1])):
                    cur = _ddb_cx(b, cur, c, t)
            else:
                raise ValueError(f"dd witness: unsupported gate {name!r}")
        else:
            raise ValueError(f"dd witness: unsupported gate {name!r}")
    return cur


def _ddb_identity(b, n):
    cur = (1.0 + 0.0j, 0)
    for lv in range(n):
        cur = b.mk(lv, (1.0 + 0.0j, cur[1], 0.0 + 0.0j, 0,
                        0.0 + 0.0j, 0, 1.0 + 0.0j, cur[1]))
    return cur


def _ddb_cx(b, a, c, t):
    h = (complex(0.0, 0.7071067811865476), complex(0.0, 0.7071067811865476),
         complex(0.0, 0.7071067811865476), complex(0.0, -0.7071067811865476))
    a = b.apply_1q(a, t, h)
    a = b.apply_cz(a, max(c, t), min(c, t), -1.0 + 0.0j)
    a = b.apply_1q(a, t, h)
    return a


def dd_equal(ops_a, ops_b, n, tol=1e-7, max_nodes=400_000):
    """Exact equivalence up to global phase via independent decision
    diagrams.  Returns None when the diagram exceeds the checker's own
    node budget (the checker never guesses)."""
    b = _DDB()
    try:
        ra = _dd_ops_unitary(b, ops_a, n)
        rb = _dd_ops_unitary(b, ops_b, n)
    except RecursionError:
        return None
    if len(b.nodes) > max_nodes:
        return None
    val = b.overlap(ra, rb, {})
    d = 1 << n
    return abs(val) / d >= 1.0 - tol
