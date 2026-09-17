"""Native two-qubit gate synthesis: emit circuits in a machine's own 2q gate.

Supported targets are the supercontrolled family U_d(pi/4, b, 0), which
contains every common native entangling gate:

    cx    (b = 0)      cz    (b = 0)      ecr   (b = 0)      iswap (b = pi/4)

The decomposer is ported from qiskit-terra 0.46 ``TwoQubitBasisDecomposer``
(Cross et al., arXiv:1811.12926), re-expressed on this repo's pure-stdlib
matrix helpers.  For each target unitary it picks the minimal gate count
(0-3) by the exact trace-fidelity table and returns the surrounding local
gates; the caller re-synthesizes locals into 1q runs.
"""
from __future__ import annotations

import cmath
import math

from .circuit import Circuit, Gate
from .kak import (mmul2, mmul4, dagger2, rz_mat, weyl, split_product,
                  _flatten2, synth_2q_ops)

# ---------------------------------------------------------------------------
# native gate unitaries (verified against qiskit Operator to 1e-15)

NATIVE_MATRICES = {
    "cx": [
        [complex(1, 0), complex(0, 0), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(1, 0), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(0, 0), complex(0, 0), complex(1, 0)],
        [complex(0, 0), complex(0, 0), complex(1, 0), complex(0, 0)],
    ],
    "cz": [
        [complex(1, 0), complex(0, 0), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(1, 0), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(0, 0), complex(1, 0), complex(0, 0)],
        [complex(0, 0), complex(0, 0), complex(0, 0), complex(-1, 0)],
    ],
    "ecr": [
        [complex(0, 0), complex(0.7071067812, 0), complex(0, 0), complex(0, 0.7071067812)],
        [complex(0.7071067812, 0), complex(0, 0), complex(0, -0.7071067812), complex(0, 0)],
        [complex(0, 0), complex(0, 0.7071067812), complex(0, 0), complex(0.7071067812, 0)],
        [complex(0, -0.7071067812), complex(0, 0), complex(0.7071067812, 0), complex(0, 0)],
    ],
    "iswap": [
        [complex(1, 0), complex(0, 0), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(0, 0), complex(0, 1), complex(0, 0)],
        [complex(0, 0), complex(0, 1), complex(0, 0), complex(0, 0)],
        [complex(0, 0), complex(0, 0), complex(0, 0), complex(1, 0)],
    ],
}
NATIVE_GATES = ("cx", "cz", "ecr", "iswap")

_ipz = [[complex(0, 1), complex(0, 0)], [complex(0, 0), complex(0, -1)]]


def _mul(a, b):
    return mmul2(a, b)


class BasisDecomposer:
    """Minimal-count synthesis of arbitrary SU(4) into a supercontrolled
    basis gate U_d(pi/4, b, 0) plus local gates (terra 0.46 port)."""

    def __init__(self, native: str):
        if native not in NATIVE_MATRICES:
            raise ValueError(f"unknown native gate {native!r}; choose from {NATIVE_GATES}")
        self.native = native
        a, b, c, K1l, K1r, K2l, K2r = weyl(NATIVE_MATRICES[native])
        if not (math.isclose(a, math.pi / 4, abs_tol=1e-9)
                and abs(c) < 1e-9):
            raise ValueError(f"native gate {native!r} is not supercontrolled")
        self.basis_b = b
        self.basis_locals = (K1l, K1r, K2l, K2r)

        # Create some useful matrices U1, U2, U3 equivalent to the basis,
        # expanded as Ui = Ki1 . Ubasis . Ki2  (terra constants, function of b)
        K11l = [
            [complex(0, -1) * cmath.exp(-1j * b) / (1 + 1j), cmath.exp(-1j * b) / (1 + 1j)],
            [complex(0, -1) * cmath.exp(1j * b) / (1 + 1j), -cmath.exp(1j * b) / (1 + 1j)],
        ]
        K11r = [
            [complex(0, 1) * cmath.exp(-1j * b) / math.sqrt(2), -cmath.exp(-1j * b) / math.sqrt(2)],
            [cmath.exp(1j * b) / math.sqrt(2), complex(0, -1) * cmath.exp(1j * b) / math.sqrt(2)],
        ]
        K12l = [
            [complex(0, 1) / (1 + 1j), complex(0, 1) / (1 + 1j)],
            [complex(-1, 0) / (1 + 1j), complex(1, 0) / (1 + 1j)],
        ]
        K12r = [
            [complex(0, 1) / math.sqrt(2), complex(1, 0) / math.sqrt(2)],
            [complex(-1, 0) / math.sqrt(2), complex(0, -1) / math.sqrt(2)],
        ]
        K32lK21l = [
            [(1 + 1j * math.cos(2 * b)) / math.sqrt(2), 1j * math.sin(2 * b) / math.sqrt(2)],
            [1j * math.sin(2 * b) / math.sqrt(2), (1 - 1j * math.cos(2 * b)) / math.sqrt(2)],
        ]
        K21r = [
            [complex(0, -1) * cmath.exp(-2j * b) / (1 - 1j), cmath.exp(-2j * b) / (1 - 1j)],
            [complex(0, 1) * cmath.exp(2j * b) / (1 - 1j), cmath.exp(2j * b) / (1 - 1j)],
        ]
        K22l = [[complex(1 / math.sqrt(2), 0), complex(-1 / math.sqrt(2), 0)],
                [complex(1 / math.sqrt(2), 0), complex(1 / math.sqrt(2), 0)]]
        K22r = [[complex(0, 0), complex(1, 0)], [complex(-1, 0), complex(0, 0)]]
        K31l = [
            [cmath.exp(-1j * b) / math.sqrt(2), cmath.exp(-1j * b) / math.sqrt(2)],
            [-cmath.exp(1j * b) / math.sqrt(2), cmath.exp(1j * b) / math.sqrt(2)],
        ]
        K31r = [[cmath.exp(1j * b) * 1j, complex(0, 0)], [complex(0, 0), -cmath.exp(-1j * b) * 1j]]
        K32r = [
            [cmath.exp(1j * b) / (1 - 1j), -cmath.exp(-1j * b) / (1 - 1j)],
            [-1j * cmath.exp(1j * b) / (1 - 1j), -1j * cmath.exp(-1j * b) / (1 - 1j)],
        ]
        k1ld, k1rd = dagger2(K1l), dagger2(K1r)
        k2ld, k2rd = dagger2(K2l), dagger2(K2r)

        # fixed parts of the 3-gate decomposition
        self.u0l = _mul(K31l, k1ld)
        self.u0r = _mul(K31r, k1rd)
        self.u1l = _mul(_mul(k2ld, K32lK21l), k1ld)
        self.u1ra = _mul(k2rd, K32r)
        self.u1rb = _mul(K21r, k1rd)
        self.u2la = _mul(k2ld, K22l)
        self.u2lb = _mul(K11l, k1ld)
        self.u2ra = _mul(k2rd, K22r)
        self.u2rb = _mul(K11r, k1rd)
        self.u3l = _mul(k2ld, K12l)
        self.u3r = _mul(k2rd, K12r)

        # fixed parts of the 2-gate decomposition
        self.q0l = _mul(dagger2(K12l), k1ld)
        self.q0r = _mul(_mul(dagger2(K12r), _ipz), k1rd)
        self.q1la = _mul(k2ld, dagger2(K11l))
        self.q1lb = _mul(K11l, k1ld)
        self.q1ra = _mul(_mul(k2rd, _ipz), dagger2(K11r))
        self.q1rb = _mul(K11r, k1rd)
        self.q2l = _mul(k2ld, K12l)
        self.q2r = _mul(k2rd, K12r)

    # -- fidelity table ----------------------------------------------------
    def traces(self, ta, tb, tc):
        """|Tr| of the best 0/1/2/3-gate approximations of U_d(ta, tb, tc)."""
        bb = self.basis_b
        p4 = math.pi / 4
        return [
            4 * abs(complex(math.cos(ta) * math.cos(tb) * math.cos(tc),
                            math.sin(ta) * math.sin(tb) * math.sin(tc))),
            4 * abs(complex(math.cos(p4 - ta) * math.cos(bb - tb) * math.cos(tc),
                            math.sin(p4 - ta) * math.sin(bb - tb) * math.sin(tc))),
            4 * abs(math.cos(tc)),
            4.0,
        ]

    # -- decompositions ----------------------------------------------------
    @staticmethod
    def _decomp0(target):
        K1l, K1r, K2l, K2r = target
        return _mul(K1l, K2l), _mul(K1r, K2r)

    def _decomp1(self, target, basis):
        bK1l, bK1r, bK2l, bK2r = basis
        K1l, K1r, K2l, K2r = target
        u0l, u0r = _mul(K1l, dagger2(bK1l)), _mul(K1r, dagger2(bK1r))
        u1l, u1r = _mul(dagger2(bK2l), K2l), _mul(dagger2(bK2r), K2r)
        return u1r, u1l, u0r, u0l

    def _decomp2(self, target, rz):
        K1l, K1r, K2l, K2r = target
        ta, tb = self._ta, self._tb
        u0l, u0r = _mul(K1l, self.q0l), _mul(K1r, self.q0r)
        u1l = _mul(_mul(self.q1la, rz(-2 * ta)), self.q1lb)
        u1r = _mul(_mul(self.q1ra, rz(2 * tb)), self.q1rb)
        u2l, u2r = _mul(self.q2l, K2l), _mul(self.q2r, K2r)
        return u2r, u2l, u1r, u1l, u0r, u0l

    def _decomp3(self, target, rz):
        K1l, K1r, K2l, K2r = target
        ta, tb, tc = self._ta, self._tb, self._tc
        u0l, u0r = _mul(K1l, self.u0l), _mul(K1r, self.u0r)
        u1l = self.u1l
        u1r = _mul(_mul(self.u1ra, rz(-2 * tc)), self.u1rb)
        u2l = _mul(_mul(self.u2la, rz(-2 * ta)), self.u2lb)
        u2r = _mul(_mul(self.u2ra, rz(2 * tb)), self.u2rb)
        u3l, u3r = _mul(self.u3l, K2l), _mul(self.u3r, K2r)
        return u3r, u3l, u2r, u2l, u1r, u1l, u0r, u0l

    # -- public ------------------------------------------------------------
    def synth(self, U4, fid_tol: float = 1 - 1e-9):
        """Return (n_basis_gates, [local 2x2 layers]) with n minimal such that
        the fidelity of the synthesis is >= fid_tol.

        Local layers come back as [(l0, r0), (l1, r1), ...] in circuit order
        (the native gate sits between consecutive layers).
        """
        ta, tb, tc, K1l, K1r, K2l, K2r = weyl(U4)
        self._ta, self._tb, self._tc = ta, tb, tc
        target = (K1l, K1r, K2l, K2r)
        basis = self.basis_locals
        tr = self.traces(ta, tb, tc)
        n = next((k for k in range(4) if tr[k] / 4.0 >= fid_tol), 3)
        if n == 0:
            ll, rr = self._decomp0(target)
            return 0, [(ll, rr)]
        if n == 1:
            u1r, u1l, u0r, u0l = self._decomp1(target, basis)
            return 1, [(u0l, u0r), (u1l, u1r)]
        if n == 2:
            u2r, u2l, u1r, u1l, u0r, u0l = self._decomp2(target, rz_mat)
            return 2, [(u0l, u0r), (u1l, u1r), (u2l, u2r)]
        u3r, u3l, u2r, u2l, u1r, u1l, u0r, u0l = self._decomp3(target, rz_mat)
        return 3, [(u0l, u0r), (u1l, u1r), (u2l, u2r), (u3l, u3r)]


_DECOMP: dict = {}


def _decomposer(native: str) -> BasisDecomposer:
    d = _DECOMP.get(native)
    if d is None:
        d = _DECOMP[native] = BasisDecomposer(native)
    return d


def _locals_to_ops(loc, rloc, qa: int, qb: int, prefer_u: bool = True):
    # convention note (empirically pinned against terra): our weyl locals
    # place `loc` on the high-bit wire (qb) and `rloc` on the low wire (qa)
    from .linalg import resynth
    ops = []
    lb = resynth(_flatten2(loc), prefer_u=prefer_u)
    ops += [Gate(g.name, g.params, (qb,)) for g in lb]
    la = resynth(_flatten2(rloc), prefer_u=prefer_u)
    ops += [Gate(g.name, g.params, (qa,)) for g in la]
    return ops


def synth_2q_native(U4, native: str, fid_tol: float = 1 - 1e-9,
                    prefer_u: bool = True) -> list:
    """Synthesize a 4x4 unitary into {native gate + 1q locals}.

    Returned gates act on wires (0, 1) for the native gate; locals carry
    wire indices 0/1 already.  The caller remaps to real qubit pairs.
    """
    from .linalg import resynth
    d = _decomposer(native)
    n, layers = d.synth(U4, fid_tol=fid_tol)
    ops = []
    for i, (loc, rloc) in enumerate(reversed(layers)):
        ops += _locals_to_ops(loc, rloc, 0, 1, prefer_u)
        if i < n:
            ops.append(Gate(native, (), (0, 1)))
    return ops


def rebase(circ: Circuit, native: str, oneq: str = "u3",
           fid_tol: float = 1 - 1e-9) -> Circuit:
    """Re-express every two-qubit block of `circ` in the native gate basis.

    `oneq`: "u3" keeps single-qubit runs as u3; "rz-sx-x" translates them to
    the IBM hardware 1q basis.  Gates outside blocks pass through; loose
    cx/cz/cp gates are treated as 1-gate blocks and synthesized too.
    """
    from .kak import collect_blocks, _block_circuit
    if native not in NATIVE_MATRICES:
        raise ValueError(f"unknown native gate {native!r}; choose from {NATIVE_GATES}")
    d = _decomposer(native)
    ops = list(circ.ops)
    blocks = collect_blocks(circ)
    out = []
    prev = 0
    for (start, end, a, b) in blocks:
        bc = _block_circuit(circ, start, end, a, b)
        U4 = None
        from .equivalence import unitary as _unitary
        try:
            U4 = _unitary(bc)
        except Exception:
            U4 = None
        if U4 is None:
            continue
        n, layers = d.synth(U4, fid_tol=fid_tol)
        out.extend(ops[prev:start])
        for i, (loc, rloc) in enumerate(reversed(layers)):
            out += _locals_to_ops(loc, rloc, a, b)
            if i < n:
                out.append(Gate(native, (), (a, b)))
        prev = end
    out.extend(ops[prev:])
    res = Circuit(circ.num_qubits, out)
    if oneq == "rz-sx-x":
        from .hardware import translate_1q_to_rz_sx_x
        res = translate_1q_to_rz_sx_x(res)
        res = _squash_rzsxx(res)
    elif oneq == "u3":
        from .optimize import u3_fold
        res = u3_fold(res)
    return res


def _squash_rzsxx(circ: Circuit) -> Circuit:
    """Exact 1q cleanup for the {rz, sx, x} basis: merge adjacent rz runs,
    sx.sx -> x, x.x -> id (iterated to fixpoint)."""
    changed = True
    ops = list(circ.ops)
    while changed:
        changed = False
        out = []
        i = 0
        while i < len(ops):
            g = ops[i]
            if i + 1 < len(ops):
                h = ops[i + 1]
                same_w = g.qubits == h.qubits and len(g.qubits) == 1
                if same_w and g.name == "rz" and h.name == "rz":
                    a = g.params[0] if g.params else 0.0
                    b2 = h.params[0] if h.params else 0.0
                    ang = (a + b2 + math.pi) % (2 * math.pi) - math.pi
                    if abs(ang) > 1e-12:
                        out.append(Gate("rz", (ang,), g.qubits))
                    i += 2
                    changed = True
                    continue
                if same_w and g.name == "sx" and h.name == "sx":
                    out.append(Gate("x", (), g.qubits))
                    i += 2
                    changed = True
                    continue
                if same_w and g.name == "x" and h.name == "x":
                    i += 2
                    changed = True
                    continue
            out.append(g)
            i += 1
        ops = out
    return Circuit(circ.num_qubits, ops)


def _native_expansion(native: str):
    """Decompose a native gate into {u3 + cx} ops on wires (0, 1)."""
    return synth_2q_ops(NATIVE_MATRICES[native], fid_tol=1 - 1e-12)


NATIVE_EXPANSION = {name: _native_expansion(name) for name in ("cz", "ecr", "iswap")}
