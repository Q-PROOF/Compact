"""Two-qubit KAK / Weyl decomposition + minimal-CX re-synthesis (pure Python).

Numerically-stable algorithm from Cross et al., arXiv:1811.12926 App. B (the
same algorithm Qiskit's TwoQubitWeylDecomposition used before its Rust port):
simultaneous diagonalization of the real and imaginary parts of the
complex-symmetric magic-basis matrix M2 = Up^T Up, followed by the
Weyl-chamber canonicalization flips (Kraus & Cirac 2001).

Everything here is exact up to global phase; global phase is dropped because
every consumer (fidelity checks, benchmarks) is phase-insensitive.

Safety contract: every re-synthesized block is verified against the block's
own 4x4 unitary (up to global phase) before replacing it.  The caller
compares the merged result against the input and keeps the best.
"""
from __future__ import annotations

import cmath
import math

from .circuit import Circuit, Gate
from .equivalence import check_equivalent
from .linalg import resynth

try:  # optional Rust acceleration; the pure-Python path is always the fallback
    import compactq_native as _NATIVE
except Exception:
    _NATIVE = None

_PI = math.pi
_PI2 = _PI / 2
_PI4 = _PI / 4

# ---------------------------------------------------------------------------
# 4x4 complex linear algebra (rows = list of lists of complex)
# ---------------------------------------------------------------------------

def mmul4(A, B):
    out = [[0j] * 4 for _ in range(4)]
    for i in range(4):
        Ai = A[i]
        Oi = out[i]
        for k in range(4):
            aik = Ai[k]
            if aik == 0:
                continue
            Bk = B[k]
            Oi[0] += aik * Bk[0]
            Oi[1] += aik * Bk[1]
            Oi[2] += aik * Bk[2]
            Oi[3] += aik * Bk[3]
    return out


def kron22(L, R):
    """Kronecker product L (x) R; L acts on the HIGH bit (qubit 1)."""
    out = [[0j] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            lij = L[i][j]
            if lij == 0:
                continue
            for r in range(2):
                for c in range(2):
                    out[2 * i + r][2 * j + c] = lij * R[r][c]
    return out


def transpose4(A):
    return [[A[j][i] for j in range(4)] for i in range(4)]


def det4(A):
    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    total = 0j
    for j in range(4):
        minor = [[A[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        term = A[0][j] * det3(minor)
        total += term if j % 2 == 0 else -term
    return total


def det_real4(A):
    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    total = 0.0
    for j in range(4):
        minor = [[A[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        term = A[0][j] * det3(minor)
        total += term if j % 2 == 0 else -term
    return total


# 2x2 algebra (rows = list of lists of complex)
def mmul2(A, B):
    return [[A[i][0] * B[0][j] + A[i][1] * B[1][j] for j in range(2)]
            for i in range(2)]


def dagger2(A):
    return [[A[j][i].conjugate() for j in range(2)] for i in range(2)]


def rz_mat(theta):
    e = cmath.exp(-0.5j * theta)
    return [[e, 0j], [0j, e.conjugate()]]


_ipx = [[0j, 1j], [1j, 0j]]
_ipy = [[0j, 1.0], [-1.0, 0j]]
_ipz = [[1j, 0j], [0j, -1j]]

# Magic basis (unnormalized, as in qiskit-terra: B @ B^dagger = 4*I/... the
# 0.5 factor on the dagger makes the pair exactly unitary).
_B = [[1.0, 1j, 0j, 0j],
      [0j, 0j, 1j, 1.0],
      [0j, 0j, 1j, -1.0],
      [1.0, -1j, 0j, 0j]]
_Bd = [[0.5 * _B[j][i].conjugate() for j in range(4)] for i in range(4)]


def magic_fwd(U):
    return mmul4(mmul4(_B, U), _Bd)


def magic_rev(U):
    return mmul4(mmul4(_Bd, U), _B)


# ---------------------------------------------------------------------------
# Jacobi eigenvalue algorithm for real symmetric 4x4 matrices
# ---------------------------------------------------------------------------

def jacobi_eigh(A):
    """Diagonalize real symmetric A; returns (eigenvalues, V) with A ~= V D V^T."""
    n = 4
    a = [row[:] for row in A]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(50):
        off = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                off += a[i][j] * a[i][j]
        if off < 1e-30:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                theta = (a[q][q] - a[p][p]) / (2.0 * apq)
                sign = 1.0 if theta >= 0 else -1.0
                t = sign / (abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    return [a[i][i] for i in range(n)], v


# deterministic pseudo-random mix pairs for degenerate M2 retries
def _mix_pairs(count=12):
    state = 2020
    pairs = []
    for _ in range(count):
        state = (6364136223846793005 * state + 1442695040888963407) % (1 << 64)
        w1 = ((state >> 11) / float(1 << 53)) * 2.0 - 1.0
        state = (6364136223846793005 * state + 1442695040888963407) % (1 << 64)
        w2 = ((state >> 11) / float(1 << 53)) * 2.0 - 1.0
        pairs.append((w1, w2))
    return pairs


# ---------------------------------------------------------------------------
# SU(2) x SU(2) product-gate split
# ---------------------------------------------------------------------------

def split_product(M):
    """M = L (x) R with L, R in SU(2); returns (L, R) or None."""
    R2 = [[M[0][0], M[0][1]], [M[1][0], M[1][1]]]
    detR = R2[0][0] * R2[1][1] - R2[0][1] * R2[1][0]
    if abs(detR) < 0.1:
        R2 = [[M[2][0], M[2][1]], [M[3][0], M[3][1]]]
        detR = R2[0][0] * R2[1][1] - R2[0][1] * R2[1][0]
    if abs(detR) < 0.1:
        return None
    s = cmath.sqrt(detR)
    R2 = [[x / s for x in row] for row in R2]

    # temp = M @ (I (x) R^dagger); strided block = L
    Rd = dagger2(R2)
    ir = [[Rd[0][0], Rd[0][1], 0j, 0j],
          [Rd[1][0], Rd[1][1], 0j, 0j],
          [0j, 0j, Rd[0][0], Rd[0][1]],
          [0j, 0j, Rd[1][0], Rd[1][1]]]
    temp = mmul4(M, ir)
    L2 = [[temp[0][0], temp[0][2]], [temp[2][0], temp[2][2]]]
    detL = L2[0][0] * L2[1][1] - L2[0][1] * L2[1][0]
    if abs(detL) < 0.9:
        return None
    s = cmath.sqrt(detL)
    L2 = [[x / s for x in row] for row in L2]
    return L2, R2


# ---------------------------------------------------------------------------
# Weyl decomposition
# ---------------------------------------------------------------------------

def weyl(U4):
    """Weyl decomposition of a 4x4 unitary.

    Returns (a, b, c, K1l, K1r, K2l, K2r) with
    U = (K1l x K1r) . exp(i(a XX + b YY + c ZZ)) . (K2l x K2r)
    up to global phase, pi/4 >= a >= b >= |c|;  or None on failure.
    """
    detU = det4(U4)
    if abs(detU) < 1e-12:
        return None
    root = detU ** 0.25
    U = [[x / root for x in row] for row in U4]

    Up = magic_rev(U)
    M2 = mmul4(transpose4(Up), Up)

    P = None
    D = None
    for w1, w2 in _mix_pairs():
        A = [[(w1 * M2[i][j].real + w2 * M2[i][j].imag) for j in range(4)]
             for i in range(4)]
        _, Ptry = jacobi_eigh(A)
        # D_i = (P^T M2 P)[i][i] = sum_kl P[k][i] * M2[k][l] * P[l][i]
        D = [sum(Ptry[k][i] * M2[k][l] * Ptry[l][i]
                 for k in range(4) for l in range(4))
             for i in range(4)]
        # reconstruction check: P diag(D) P^T ~= M2
        PD = [[Ptry[i][j] * D[j] for j in range(4)] for i in range(4)]
        recon = mmul4(PD, transpose4(Ptry))
        scale = max(1.0, max(abs(M2[i][j]) for i in range(4) for j in range(4)))
        ok = all(abs(recon[i][j] - M2[i][j]) < 1e-11 * scale
                 for i in range(4) for j in range(4))
        if ok:
            P = Ptry
            break
    if P is None:
        return None

    d = [-cmath.phase(x) / 2 for x in D]
    d[3] = -d[0] - d[1] - d[2]
    cs = [((d[k] + d[3]) / 2) % (2 * _PI) for k in range(3)]

    # Reorder the eigenvalues to get in the Weyl chamber
    cstemp = [x % _PI2 for x in cs]
    cstemp = [min(x, _PI2 - x) for x in cstemp]
    order = sorted(range(3), key=cstemp.__getitem__)
    order = [order[1], order[2], order[0]]
    cs = [cs[k] for k in order]
    d[:3] = [d[k] for k in order]
    P = [[P[i][order[j]] if j < 3 else P[i][3] for j in range(4)] for i in range(4)]

    # Fix the sign of P to be in SO(4)
    if det_real4(P) < 0:
        P = [row[:] for row in P]
        for i in range(4):
            P[i][3] = -P[i][3]

    E = [[cmath.exp(1j * d[j]) if i == j else 0j for j in range(4)] for i in range(4)]
    K1 = magic_fwd(mmul4(mmul4(Up, P), E))
    K2 = magic_fwd(transpose4(P))

    sp1 = split_product(K1)
    sp2 = split_product(K2)
    if sp1 is None or sp2 is None:
        return None
    K1l, K1r = sp1
    K2l, K2r = sp2

    # Flip into Weyl chamber
    if cs[0] > _PI2:
        cs[0] -= 3 * _PI2
        K1l = mmul2(K1l, _ipy)
        K1r = mmul2(K1r, _ipy)
    if cs[1] > _PI2:
        cs[1] -= 3 * _PI2
        K1l = mmul2(K1l, _ipx)
        K1r = mmul2(K1r, _ipx)
    conjs = 0
    if cs[0] > _PI4:
        cs[0] = _PI2 - cs[0]
        K1l = mmul2(K1l, _ipy)
        K2r = mmul2(_ipy, K2r)
        conjs += 1
    if cs[1] > _PI4:
        cs[1] = _PI2 - cs[1]
        K1l = mmul2(K1l, _ipx)
        K2r = mmul2(_ipx, K2r)
        conjs += 1
    if cs[2] > _PI2:
        cs[2] -= 3 * _PI2
        K1l = mmul2(K1l, _ipz)
        K1r = mmul2(K1r, _ipz)
    if conjs == 1:
        cs[2] = _PI2 - cs[2]
        K1l = mmul2(K1l, _ipz)
        K2r = mmul2(_ipz, K2r)
    if cs[2] > _PI4:
        cs[2] -= _PI2
        K1l = mmul2(K1l, _ipz)
        K1r = mmul2(K1r, _ipz)

    return cs[1], cs[0], cs[2], K1l, K1r, K2l, K2r


# ---------------------------------------------------------------------------
# CX-basis circuit templates (ported from qiskit-terra TwoQubitBasisDecomposer)
# ---------------------------------------------------------------------------

def _cx_matrix():
    return [[1.0, 0j, 0j, 0j],
            [0j, 0j, 0j, 1.0],
            [0j, 0j, 1.0, 0j],
            [0j, 1.0, 0j, 0j]]


_BASIS = None  # cached weyl(CX)


def _basis():
    global _BASIS
    if _BASIS is None:
        res = weyl(_cx_matrix())
        if res is None:
            raise RuntimeError("KAK bootstrap failed on the CX gate")
        _BASIS = res
    return _BASIS


def _decomposer_mats():
    """Fixed matrices for the CX-basis decomp2/decomp3 templates."""
    _, _, basis_b, bK1l, bK1r, bK2l, bK2r = _basis()
    b = basis_b
    k1ld, k1rd = dagger2(bK1l), dagger2(bK1r)
    k2ld, k2rd = dagger2(bK2l), dagger2(bK2r)

    def c(exp_i_theta):
        return cmath.exp(1j * exp_i_theta)

    K11l = [[c(-b) * -1j / (1 + 1j), c(-b) / (1 + 1j)],
            [c(b) * -1j / (1 + 1j), -c(b) / (1 + 1j)]]
    K11r = [[1j * c(-b) / math.sqrt(2), -c(-b) / math.sqrt(2)],
            [c(b) / math.sqrt(2), -1j * c(b) / math.sqrt(2)]]
    K12l = [[1j / (1 + 1j), 1j / (1 + 1j)], [-1 / (1 + 1j), 1 / (1 + 1j)]]
    K12r = [[1j / math.sqrt(2), 1 / math.sqrt(2)],
            [-1 / math.sqrt(2), -1j / math.sqrt(2)]]
    cb, sb = math.cos(2 * b), math.sin(2 * b)
    K32lK21l = [[(1 + 1j * cb) / math.sqrt(2), 1j * sb / math.sqrt(2)],
                [1j * sb / math.sqrt(2), (1 - 1j * cb) / math.sqrt(2)]]
    K21r = [[c(-2 * b) * -1j / (1 - 1j), c(-2 * b) / (1 - 1j)],
            [c(2 * b) * 1j / (1 - 1j), c(2 * b) / (1 - 1j)]]
    K22l = [[1 / math.sqrt(2), -1 / math.sqrt(2)],
            [1 / math.sqrt(2), 1 / math.sqrt(2)]]
    K22r = [[0j, 1.0], [-1.0, 0j]]
    K31l = [[c(-b) / math.sqrt(2), c(-b) / math.sqrt(2)],
            [-c(b) / math.sqrt(2), c(b) / math.sqrt(2)]]
    K31r = [[1j * c(b), 0j], [0j, -1j * c(-b)]]
    K32r = [[c(b) / (1 - 1j), -c(-b) / (1 - 1j)],
            [-1j * c(b) / (1 - 1j), -1j * c(-b) / (1 - 1j)]]

    m = {}
    m["u0l"] = mmul2(K31l, k1ld)
    m["u0r"] = mmul2(K31r, k1rd)
    m["u1l"] = mmul2(mmul2(k2ld, K32lK21l), k1ld)
    m["u1ra"] = mmul2(k2rd, K32r)
    m["u1rb"] = mmul2(K21r, k1rd)
    m["u2la"] = mmul2(k2ld, K22l)
    m["u2lb"] = mmul2(K11l, k1ld)
    m["u2ra"] = mmul2(k2rd, K22r)
    m["u2rb"] = mmul2(K11r, k1rd)
    m["u3l"] = mmul2(k2ld, K12l)
    m["u3r"] = mmul2(k2rd, K12r)

    m["q0l"] = mmul2(dagger2(K12l), k1ld)
    m["q0r"] = mmul2(mmul2(dagger2(K12r), _ipz), k1rd)
    m["q1la"] = mmul2(k2ld, dagger2(K11l))
    m["q1lb"] = mmul2(K11l, k1ld)
    m["q1ra"] = mmul2(mmul2(k2rd, _ipz), dagger2(K11r))
    m["q1rb"] = mmul2(K11r, k1rd)
    m["q2l"] = mmul2(k2ld, K12l)
    m["q2r"] = mmul2(k2rd, K12r)
    return m


_MATS = None


def _mats():
    global _MATS
    if _MATS is None:
        _MATS = _decomposer_mats()
    return _MATS


def _fid_table(a, b, c):
    """Average-gate fidelity of the best nbasis-CX approximation, k = 0..3."""
    ta, tb, tc = a, b, c
    t0 = 4 * complex(math.cos(ta) * math.cos(tb) * math.cos(tc),
                     math.sin(ta) * math.sin(tb) * math.sin(tc))
    t1 = 4 * complex(math.cos(_PI4 - ta) * math.cos(-tb) * math.cos(tc),
                     math.sin(_PI4 - ta) * math.sin(-tb) * math.sin(tc))
    t2 = 4 * math.cos(tc)
    t3 = 4.0 + 0j
    return [abs(t0) / 4, abs(t1) / 4, abs(t2) / 4, abs(t3) / 4]


def _decomp_locals(nb, a, b, c, K1l, K1r, K2l, K2r):
    """Local 2x2 layers in circuit-time order interleaved with CX markers.

    Returns a list of items: ('l', matrix, wire) or ('cx',).
    wire 0 = low qubit, wire 1 = high qubit (little-endian, matches unitary()).
    """
    m = _mats()
    rz = rz_mat
    if nb == 3:
        U0l = mmul2(K1l, m["u0l"])
        U0r = mmul2(K1r, m["u0r"])
        U1l = m["u1l"]
        U1r = mmul2(mmul2(m["u1ra"], rz(-2 * c)), m["u1rb"])
        U2l = mmul2(mmul2(m["u2la"], rz(-2 * a)), m["u2lb"])
        U2r = mmul2(mmul2(m["u2ra"], rz(2 * b)), m["u2rb"])
        U3l = mmul2(m["u3l"], K2l)
        U3r = mmul2(m["u3r"], K2r)
        return [("l", U3r, 0), ("l", U3l, 1), ("cx",),
                ("l", U2r, 0), ("l", U2l, 1), ("cx",),
                ("l", U1r, 0), ("l", U1l, 1), ("cx",),
                ("l", U0r, 0), ("l", U0l, 1)]
    if nb == 2:
        U0l = mmul2(K1l, m["q0l"])
        U0r = mmul2(K1r, m["q0r"])
        U1l = mmul2(mmul2(m["q1la"], rz(-2 * a)), m["q1lb"])
        U1r = mmul2(mmul2(m["q1ra"], rz(2 * b)), m["q1rb"])
        U2l = mmul2(m["q2l"], K2l)
        U2r = mmul2(m["q2r"], K2r)
        return [("l", U2r, 0), ("l", U2l, 1), ("cx",),
                ("l", U1r, 0), ("l", U1l, 1), ("cx",),
                ("l", U0r, 0), ("l", U0l, 1)]
    if nb == 1:
        bK1l, bK1r, bK2l, bK2r = _basis()[3:]
        U0l = mmul2(K1l, dagger2(bK1l))
        U0r = mmul2(K1r, dagger2(bK1r))
        U1l = mmul2(dagger2(bK2l), K2l)
        U1r = mmul2(dagger2(bK2r), K2r)
        return [("l", U1r, 0), ("l", U1l, 1), ("cx",),
                ("l", U0r, 0), ("l", U0l, 1)]
    # nb == 0: pure product gate
    U0l = mmul2(K1l, K2l)
    U0r = mmul2(K1r, K2r)
    return [("l", U0r, 0), ("l", U0l, 1)]


def _flatten2(L):
    return (L[0][0], L[0][1], L[1][0], L[1][1])


def synth_2q_ops(U4, fid_tol=1 - 1e-13, force_nb=None):
    """Synthesize a 4x4 unitary into gates on wires (0, 1) with minimal CX.

    Returns a list of Gate objects or None when synthesis fails.
    The result is exact up to global phase (fid_tol allows <1e-13 slack).
    force_nb pins the basis-gate count (approximate when < 3) instead of
    choosing by the fidelity table.
    """
    try:
        w = weyl(U4)
        if w is None:
            return None
        a, b, c, K1l, K1r, K2l, K2r = w
        fids = _fid_table(a, b, c)
        if force_nb is not None:
            nb = max(0, min(3, int(force_nb)))
        else:
            nb = 3
            for k in range(4):
                if fids[k] >= fid_tol:
                    nb = k
                    break
        layers = _decomp_locals(nb, a, b, c, K1l, K1r, K2l, K2r)
        ops = []
        for item in layers:
            if item[0] == "cx":
                ops.append(Gate("cx", (), (0, 1)))
            else:
                _, L, wire = item
                # NB: prefer_u=False on purpose — u3 locals would block the
                # commutation/CP passes that run after KAK in the pipeline.
                # The final u3-fold pass (search.py) harvests gate count.
                for g in resynth(_flatten2(L)):
                    ops.append(Gate(g.name, g.params, (wire,)))
        return ops
    except Exception:
        return None


# ---------------------------------------------------------------------------
# block collection (maximal 2-qubit-confined segments)
# ---------------------------------------------------------------------------

def collect_blocks(circ: Circuit):
    """Maximal segments acting on exactly one qubit pair {a, b} (a < b),
    including immediately-preceding 1-qubit gates on the pair."""
    ops = circ.ops
    n = len(ops)
    used = [False] * n
    blocks = []
    i = 0
    while i < n:
        if not used[i] and len(ops[i].qubits) == 2:
            a, b = sorted(ops[i].qubits)
            j = i
            while j < n and not used[j] and set(ops[j].qubits) <= {a, b}:
                j += 1
            start = i
            while start > 0 and (not used[start - 1]) and len(ops[start - 1].qubits) == 1 \
                    and ops[start - 1].qubits[0] in (a, b):
                start -= 1
            for k in range(start, j):
                used[k] = True
            blocks.append((start, j, a, b))
            i = j
        else:
            i += 1
    return blocks


def _block_circuit(circ: Circuit, start: int, end: int, a: int, b: int) -> Circuit:
    remap = {a: 0, b: 1}
    return Circuit(2, [Gate(g.name, g.params, tuple(remap[q] for q in g.qubits))
                       for g in circ.ops[start:end]])


def _score(c: Circuit):
    return (c.two_qubit_count(), len(c.ops), c.depth())


# --------------------------------------------------------------------------
# optional external synthesis fallbacks
# --------------------------------------------------------------------------

# A pytket extension that fails to initialize once can abort the whole
# process on retry (nanobind refuses duplicate type registration), so a
# failed import is remembered and never re-attempted.
_PYTKET_DEAD = False


def _tket_squash(ops, a: int, b: int):
    """Re-synthesize a block via TKET's TwoQubitSquash (KAK-optimal CX count)."""
    global _PYTKET_DEAD
    if _PYTKET_DEAD:
        return None
    try:
        from .qiskit_bridge import to_qiskit
        from pytket.extensions.qiskit import qiskit_to_tk, tk_to_qiskit
        from pytket.passes import TwoQubitSquash
        from pytket.predicates import CompilationUnit
    except Exception:
        _PYTKET_DEAD = True
        return None
    remap = {a: 0, b: 1}
    block = Circuit(2, [Gate(g.name, g.params, tuple(remap[q] for q in g.qubits))
                        for g in ops])
    try:
        tk_block = qiskit_to_tk(to_qiskit(block))
    except Exception:
        return None
    cu = CompilationUnit(tk_block)
    try:
        if not TwoQubitSquash(allow_swaps=True).apply(cu):
            return None
    except Exception:
        return None
    out_qc = tk_to_qiskit(cu.circuit)
    out = []
    for inst in out_qc.data:
        name = inst.operation.name
        if name in ("barrier", "measure", "reset", "id"):
            continue
        qs = tuple(out_qc.find_bit(q).index for q in inst.qubits)
        params = tuple(float(p) for p in inst.operation.params)
        if len(qs) == 1 and name not in ("h", "x", "y", "z", "s", "sdg", "t", "tdg",
                                         "rx", "ry", "rz", "p", "sx", "u", "u3"):
            return None
        out.append(Gate(name, params, qs))
    return out


def _qiskit_synth(ops, a: int, b: int):
    """Re-synthesize a block via Qiskit's CX-basis decomposer -> ops or None."""
    try:
        import numpy as np
        from qiskit.synthesis import TwoQubitBasisDecomposer
        from qiskit.circuit.library import CXGate
    except Exception:
        return None
    remap = {a: 0, b: 1}
    block = Circuit(2, [Gate(g.name, g.params, tuple(remap[q] for q in g.qubits))
                        for g in ops])
    U = _block_unitary(block)
    if U is None:
        return None
    try:
        dec = TwoQubitBasisDecomposer(CXGate())
        qc = dec(np.array(U, dtype=complex))
    except Exception:
        return None  # qiskit's own decomposer can fail on degenerate inputs
    out = []
    for inst in qc.data:
        name = inst.operation.name
        if name in ("barrier", "measure", "reset", "id"):
            continue
        qs = tuple(qc.find_bit(q).index for q in inst.qubits)
        params = tuple(float(p) for p in inst.operation.params)
        if len(qs) == 1 and name not in ("h", "x", "y", "z", "s", "sdg", "t", "tdg",
                                         "rx", "ry", "rz", "p", "sx", "u", "u3"):
            return None
        out.append(Gate(name, params, qs))
    return out


def _block_unitary(block_circ: Circuit):
    """4x4 unitary of a block circuit; uses the optional compactq_native Rust
    kernels when installed (validated bit-compatible with the Python path,
    which remains the fallback and the verification reference)."""
    if _NATIVE is not None:
        try:
            names, pf, pc, qf, qc = [], [], [], [], []
            for g in block_circ.ops:
                names.append(g.name)
                pf.extend(g.params)
                pc.append(len(g.params))
                qf.extend(g.qubits)
                qc.append(len(g.qubits))
            flat = _NATIVE.block_unitary(block_circ.num_qubits, names, pf, pc, qf, qc)
            dim = 1 << block_circ.num_qubits
            return [[complex(flat[2 * (i * dim + j)], flat[2 * (i * dim + j) + 1])
                     for j in range(dim)]
                    for i in range(dim)]
        except Exception:
            pass
    try:
        from .equivalence import unitary as _unitary
        return _unitary(block_circ)
    except Exception:
        return None



# --------------------------------------------------------------------------
# cross-pair commutative merging
# --------------------------------------------------------------------------
def merge_crosspair(circ: Circuit) -> Circuit:
    """Move gates that are disjoint from a 2-qubit block past the block, so
    later rounds on the SAME pair merge into one KAK window.

    Gates sharing a wire with the block are blockers and end the window.
    The reorder is exact (disjoint-wire gates commute); for <= 6 qubits the
    rewritten circuit is fidelity-verified against the original as a whole
    before being returned.
    """
    ops = circ.ops
    n_ops = len(ops)
    new_ops = ops[:]
    changed = False
    used = [False] * n_ops
    i = 0
    while i < n_ops:
        g = ops[i]
        if used[i] or len(g.qubits) != 2:
            i += 1
            continue
        a, b = g.qubits
        wires = {a, b}
        block = [i]
        skip = []
        j = i + 1
        while j < n_ops:
            gj = ops[j]
            gj_set = set(gj.qubits)
            if gj_set <= wires:
                block.append(j)
            elif gj_set.isdisjoint(wires):
                skip.append(j)
            else:
                break
            j += 1
        if len([k for k in block if len(ops[k].qubits) == 2]) >= 2 and skip:
            # reorder: block first, then the commuting skips
            pos = min(block + skip)
            for off, k in enumerate(block + skip):
                new_ops[pos + off] = ops[k]
                used[k] = True
            changed = True
            i = pos + len(block) + len(skip)
        else:
            for k in block:
                used[k] = True
            i = j if j > i else i + 1
    if not changed:
        return circ
    out = Circuit(circ.num_qubits, new_ops)
    if circ.num_qubits <= 6 and not check_equivalent(circ, out):
        return circ
    return out


# --------------------------------------------------------------------------
# pass
# --------------------------------------------------------------------------

def _remap_ops(ops, a: int, b: int):
    remap = {0: a, 1: b}
    return [Gate(g.name, g.params, tuple(remap[q] for q in g.qubits)) for g in ops]


def kak_pass(circ: Circuit, force: bool = False,
             fid_tol: float = 1 - 1e-13) -> Circuit:
    """Re-synthesize collected 2-qubit blocks into minimal-CX form.

    Pure-Python KAK synthesis is tried first; pytket / qiskit decomposers are
    optional fallbacks.  Every candidate is verified against the block's own
    4x4 unitary before it can replace anything.

    force=True:  emit every verified re-synthesis even if the intermediate is
                 larger (the caller merges 1q runs afterwards and compares).
    force=False: accept only blocks that are already strictly smaller.
    fid_tol:     when < 1, blocks may be replaced by a lower-CX approximation
                 whose average gate fidelity (|Tr|/d) is >= fid_tol.  The
                 exactness verification is skipped for approximate blocks.
    """
    ops = circ.ops
    blocks = collect_blocks(circ)
    if not blocks:
        return circ
    out = []
    prev = 0
    changed = False
    for (start, end, a, b) in blocks:
        out.extend(ops[prev:start])
        prev = end
        block_circ = _block_circuit(circ, start, end, a, b)
        if block_circ.two_qubit_count() < 2:
            out.extend(ops[start:end])
            continue
        U4 = _block_unitary(block_circ)
        if U4 is None:
            out.extend(ops[start:end])
            continue
        new_ops = synth_2q_ops(U4, fid_tol=fid_tol)
        if new_ops is not None:
            candidate = Circuit(2, new_ops)
        else:
            candidate = None
        if candidate is None:
            fb = _tket_squash(ops[start:end], a, b)
            if fb is None:
                fb = _qiskit_synth(ops[start:end], a, b)
            if fb is not None:
                candidate = Circuit(2, fb)
                emit = fb
            else:
                emit = None
        else:
            emit = _remap_ops(new_ops, a, b)
        if candidate is None:
            out.extend(ops[start:end])
            continue
        if fid_tol >= 1 - 1e-13:
            if not check_equivalent(block_circ, candidate):
                out.extend(ops[start:end])
                continue
        else:
            # approximate block: MEASURE the actual block fidelity (4x4,
            # cheap) and only accept when it meets the tolerance.  This makes
            # the per-block guarantee real rather than assumed from the
            # KAK class table.
            U_c = _block_unitary(candidate)
            if U_c is None:
                out.extend(ops[start:end])
                continue
                # unreachable
            tr = 0j
            for a in range(4):
                ra, rc = U4[a], U_c[a]
                for b in range(4):
                    tr += ra[b].conjugate() * rc[b]
            f = abs(tr) / 4
            if f < fid_tol:
                out.extend(ops[start:end])
                continue
        block_sc = _score(block_circ)
        cand_sc = _score(candidate)
        if force:
            # never grow the CX count; at equal CX only accept a not-larger block
            accept = (cand_sc[0] < block_sc[0]
                      or (cand_sc[0] == block_sc[0] and cand_sc[1] <= block_sc[1]))
        else:
            accept = cand_sc < block_sc
        if accept:
            out.extend(emit)
            changed = True
        else:
            out.extend(ops[start:end])
    out.extend(ops[prev:])
    if not changed:
        return circ
    return Circuit(circ.num_qubits, out)
