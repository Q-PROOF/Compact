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
    """Return G·U where G acts on `ws`.  Convention: ws[0] is the MOST
    significant local wire — matching the standard gate-matrix basis
    (|q_ws0 q_ws1 ...⟩ ordering)."""
    d = 1 << n
    shifts = [n - 1 - w for w in ws]
    m = 1 << len(ws)
    L = len(ws)
    new = [row[:] for row in U]
    for col in range(d):
        for rpat in range(m):
            acc = 0j
            for kpat in range(m):
                srow = col
                for b in range(L):
                    sh = shifts[b]
                    kb = (kpat >> (L - 1 - b)) & 1
                    srow = (srow & ~(1 << sh)) | (kb << sh)
                acc += g[rpat][kpat] * U[srow][col]
            rrow = col
            for b in range(L):
                sh = shifts[b]
                rb = (rpat >> (L - 1 - b)) & 1
                rrow = (rrow & ~(1 << sh)) | (rb << sh)
            new[rrow][col] = acc
    return new


def build_unitary(ops, n):
    d = 1 << n
    U = [[1.0 if i == j else 0.0 for j in range(d)] for i in range(d)]
    for (name, params, qubits) in ops:
        g = gate_matrix(name, params)
        U = apply_dense(U, g, qubits, n)
    return U


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
        gens.append((1 << j, 0, 0))       # X_j
    for j in range(n):
        gens.append((0, 1 << j, 0))       # Z_j
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
