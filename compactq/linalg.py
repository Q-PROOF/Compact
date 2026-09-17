"""2x2 matrix algebra and gate classification for the optimizer.

Matrices are row-major 4-tuples of complex numbers: ``(a, b, c, d)``
representing [[a, b], [c, d]].  Everything here is pure stdlib Python.
"""
from __future__ import annotations

import cmath
import math

from .circuit import Gate

try:  # optional Rust kernel: equivalent-but-faster gate decisions
    import compactq_native as _NATIVE
except Exception:
    _NATIVE = None

Mat = tuple  # (a, b, c, d) complex

I = (1 + 0j, 0j, 0j, 1 + 0j)
X = (0j, 1 + 0j, 1 + 0j, 0j)
Y = (0j, -1j, 1j, 0j)
Z = (1 + 0j, 0j, 0j, -1 + 0j)
_SQ = 1 / math.sqrt(2.0)
H = (_SQ + 0j, _SQ + 0j, _SQ + 0j, -_SQ + 0j)
S = (1 + 0j, 0j, 0j, 1j)
SDG = (1 + 0j, 0j, 0j, -1j)
T = (1 + 0j, 0j, 0j, cmath.exp(1j * math.pi / 4))
TDG = (1 + 0j, 0j, 0j, cmath.exp(-1j * math.pi / 4))

_SX = (0.5 + 0.5j, 0.5 - 0.5j, 0.5 - 0.5j, 0.5 + 0.5j)  # sqrt(X) = e^{i pi/4} RX(pi/2)
_SXDG = (0.5 - 0.5j, 0.5 + 0.5j, 0.5 + 0.5j, 0.5 - 0.5j)  # (sqrt(X))^dag
NAMED_MATRICES = {"i": I, "x": X, "y": Y, "z": Z, "h": H, "s": S, "sdg": SDG,
                  "t": T, "tdg": TDG, "sx": _SX, "sxdg": _SXDG}

# Diagonal ("phase-type") gates mapped to their equivalent P(θ) angle.
_PHASE_ANGLE = {"p": lambda p: p[0], "rz": lambda p: p[0], "s": lambda p: math.pi / 2,
                "z": lambda p: math.pi, "sdg": lambda p: -math.pi / 2,
                "t": lambda p: math.pi / 4, "tdg": lambda p: -math.pi / 4}


def mmul(A: Mat, B: Mat) -> Mat:
    a1, b1, c1, d1 = A
    a2, b2, c2, d2 = B
    return (a1 * a2 + b1 * c2, a1 * b2 + b1 * d2, c1 * a2 + d1 * c2, c1 * b2 + d1 * d2)


def mident(A: Mat) -> bool:
    a, b, c, d = A
    return (abs(a - 1) < 1e-9 and abs(b) < 1e-9 and abs(c) < 1e-9 and abs(d - 1) < 1e-9)


def same_up_to_phase(A: Mat, B: Mat, tol: float = 1e-8) -> bool:
    """True if A == e^{iφ} B for some global phase φ."""
    tr = (A[0] * B[0].conjugate() + A[1] * B[1].conjugate()
          + A[2] * B[2].conjugate() + A[3] * B[3].conjugate())
    a0, a1, a2, a3 = A
    na = (a0.real * a0.real + a0.imag * a0.imag
          + a1.real * a1.real + a1.imag * a1.imag
          + a2.real * a2.real + a2.imag * a2.imag
          + a3.real * a3.real + a3.imag * a3.imag)
    b0, b1, b2, b3 = B
    nb = (b0.real * b0.real + b0.imag * b0.imag
          + b1.real * b1.real + b1.imag * b1.imag
          + b2.real * b2.real + b2.imag * b2.imag
          + b3.real * b3.real + b3.imag * b3.imag)
    denom = na * nb
    if denom < 1e-24:
        return False
    # abs(tr)/sqrt(na*nb) > 1-tol, squared (both sides non-negative)
    return (tr.real * tr.real + tr.imag * tr.imag) > (1.0 - tol) ** 2 * denom


_GM_CACHE: dict = {}


def gate_matrix(g: Gate) -> Mat:
    """Unitary (2x2) of a 1-qubit gate (memoized on name + parameters)."""
    key = (g.name, g.params)
    hit = _GM_CACHE.get(key)
    if hit is not None:
        return hit
    m = _gate_matrix_uncached(g)
    if len(_GM_CACHE) > 8192:
        _GM_CACHE.clear()
    _GM_CACHE[key] = m
    return m


def _gate_matrix_uncached(g: Gate) -> Mat:
    n = g.name
    if n in NAMED_MATRICES:
        return NAMED_MATRICES[n]
    t = g.params[0] if g.params else 0.0
    c, s = math.cos(t / 2.0), math.sin(t / 2.0)
    if n == "rx":
        return (c + 0j, -1j * s, -1j * s, c + 0j)
    if n == "ry":
        return (c + 0j, -s + 0j, s + 0j, c + 0j)
    if n == "rz":
        e = cmath.exp(-1j * t / 2.0)
        return (e, 0j, 0j, e.conjugate())
    if n == "p":
        return (1 + 0j, 0j, 0j, cmath.exp(1j * t))
    if n in ("u", "u3"):
        tt, ph, lam = g.params
        return mat4_su2_product(rz_mat(ph), ry_mat(tt), rz_mat(lam))
    raise ValueError(f"no matrix for gate {n!r}")


def is_diagonal(g: Gate) -> bool:
    return g.name in _PHASE_ANGLE


def is_xtype(g: Gate) -> bool:
    """Gates whose conjugation by CX (as *target*) stays local and unchanged."""
    return g.name in ("x", "rx")


def wrap(x: float) -> float:
    """Wrap angle into (-pi, pi]."""
    x = math.fmod(x + math.pi, 2.0 * math.pi)
    if x <= 0:
        x += 2.0 * math.pi
    return x - math.pi


def resynth(prod: Mat, prefer_u: bool = False) -> list:
    """Resynthesise an arbitrary 1-qubit product into the fewest gates.

    Exact up to global phase.  Returns [] when the product is the identity.
    Uses the Rust `resynth_1q` kernel when available (validated equivalent
    on randomized products, identical gate counts); falls back to the
    pure-Python implementation on any error.  Strategy: named Clifford ->
    single RX/RY recognition -> exact 2-gate H-times-phase factorisations ->
    RZ-RY-RZ Euler decomposition (a single `u3` gate when `prefer_u`, used
    by the KAK synthesiser where every gate counts).
    """
    if _NATIVE is not None and hasattr(_NATIVE, "resynth_1q"):
        try:
            flat = []
            for z in prod:
                flat.extend((z.real, z.imag))
            out = _NATIVE.resynth_1q(flat, prefer_u)
            return [Gate(nm, tuple(ps), (0,)) for nm, ps in out]
        except Exception:
            pass
    return _resynth_py(prod, prefer_u)


def _resynth_py(prod: Mat, prefer_u: bool = False) -> list:
    for name, m in NAMED_MATRICES.items():
        if m is I:
            continue
        if same_up_to_phase(prod, m):
            return [Gate(name, (), (0,))]  # qubit index patched by caller

    # --- SU(2) part: V = prod / sqrt(det prod)
    a, b, c, d = prod
    det = a * d - b * c
    if abs(det) < 1e-15:
        return None  # not unitary; caller keeps original gates
    s = cmath.sqrt(det)
    va, vb, vc, vd = (x / s for x in prod)

    tol = 1e-9

    # --- single-gate RX recognition: V = [[c, -i s],[-i s, c]]
    if (abs(va - vd) < tol and abs(vb - vc) < tol
            and abs(vb.real) < tol and abs(va.imag) < tol):
        theta = 2.0 * math.atan2(-vb.imag, va.real)
        theta = wrap(theta)
        if abs(theta) < 1e-12:
            return []
        return [Gate("rx", (theta,), (0,))]

    # --- single-gate RY recognition: V real, [[c, -s],[s, c]]
    if (abs(va.imag) < tol and abs(vb.imag) < tol
            and abs(vc.imag) < tol and abs(vd.imag) < tol):
        theta = wrap(2.0 * math.atan2(vc.real, va.real))
        if abs(theta) < 1e-12:
            return []
        return [Gate("ry", (theta,), (0,))]

    # --- exact 2-gate H-times-phase factorisations.
    #  H.P(t) = (1/sqrt2)[[1, e^{it}],[1, -e^{it}]]  -> V00=V10, V01=-V11
    #  P(t).H = (1/sqrt2)[[1, 1],[e^{it}, -e^{it}]]  -> V00=V01, V10=-V11
    rsq = 1.0 / math.sqrt(2.0)
    if (abs(abs(va) - rsq) < 1e-7 and abs(abs(vb) - rsq) < 1e-7):
        if abs(va - vc) < tol and abs(vb + vd) < tol:
            phi = cmath.phase(vb) - cmath.phase(va)
            phi = wrap(phi)
            if abs(phi) < 1e-12:
                return [Gate("h", (), (0,))]
            return [Gate("p", (phi,), (0,)), Gate("h", (), (0,))]
        if abs(va - vb) < tol and abs(vc + vd) < tol:
            phi = cmath.phase(vc) - cmath.phase(va)
            phi = wrap(phi)
            if abs(phi) < 1e-12:
                return [Gate("h", (), (0,))]
            return [Gate("h", (), (0,)), Gate("p", (phi,), (0,))]

    cos2 = min(1.0, abs(va))
    theta = 2.0 * math.atan2(abs(vb), cos2)

    out: list = []
    if abs(theta) < 1e-12:
        ang = -2.0 * cmath.phase(va)  # lambda + mu
        ang = wrap(ang)
        if abs(ang) > 1e-12:
            out.append(Gate("p", (ang,), (0,)))
        return out

    # lambda - mu from the off-diagonal, lambda + mu from the diagonal.
    # NB: no 2*pi wrapping here — RZ is only 4*pi-periodic in each angle, and
    # wrapping individually would change the unitary (not just its phase).
    lam_mu = -2.0 * cmath.phase(-vb)      # lambda - mu
    if abs(cos2) > 1e-12:
        lam_mu_p = -2.0 * cmath.phase(va)  # lambda + mu
        lam = (lam_mu + lam_mu_p) / 2.0
        mu = (lam_mu_p - lam_mu) / 2.0
    else:
        lam, mu = lam_mu, 0.0

    # Derivation gives V = RZ(lam) . RY(theta) . RZ(mu).  Circuit semantics
    # apply the FIRST listed gate first (= rightmost in the matrix product),
    # so the list must be emitted in reverse: mu, theta, lam.
    if prefer_u:
        # u3(tt, ph, lam) has matrix RZ(ph) . RY(tt) . RZ(lam) -> one gate
        return [Gate("u3", (theta, lam, mu), (0,))]
    if abs(mu) > 1e-12:
        out.append(Gate("rz", (mu,), (0,)))
    out.append(Gate("ry", (theta,), (0,)))
    if abs(lam) > 1e-12:
        out.append(Gate("rz", (lam,), (0,)))
    return out


def rz_mat(t):
    e = cmath.exp(-1j * t / 2)
    return [[e, 0j], [0j, e.conjugate()]]


def ry_mat(t):
    c, s = math.cos(t / 2), math.sin(t / 2)
    return [[c, -s], [s, c]]


def mat4_su2_product(A, B, C):
    """A·B·C for 2x2 row-lists -> flattened tuple (row-major)."""
    AB = [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    M = [[sum(AB[i][k] * C[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    return (M[0][0], M[0][1], M[1][0], M[1][1])
