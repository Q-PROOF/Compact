"""Clifford stabilizer STATE simulator (CHP) - scale layer.

The density-matrix study tops out near 6 qubits and the trajectory
simulator near 14; stabilizer circuits (GHZ, BV, QAOA with Clifford
angles) run here to n ~= 60 in O(n) memory and O(gates * n) time per
shot.  This is the scale regime where closed suppression services show
their headline demos.

Rows are PYTHON BITMASKS: x[i], z[i] are ints over 2n rows (0..n-1
destabilizers, n..2n-1 stabilizers) with bit w = wire w, and r[i] the
sign bit.  Closed-form update rules for x/y/z/h/s/sdg/sx/cx/cz/swap;
any other 1q Clifford (rz at multiples of pi/2, rx(+-pi/2), ry(+-pi/2))
is conjugated NUMERICALLY - the 2x2 is built and its action on the four
Pauli types matched against signed candidates, so a wrong Y convention
cannot ship silently (this is the audit AGENTS.md asked for; the fuzz
test exercises it through random s/h/sx/rz circuits).

Noise enters as trajectory sampling: the caller injects Pauli errors per
gate (depolarizing) and readout flips per shot.  Coherent rotations and
amplitude damping are NOT representable in the Pauli formalism - the
honest scope of this layer is stochastic + readout noise at scale.
"""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate


class StabState:
    """CHP stabilizer state over `n` wires - generator rows only.

    Unlike the classic 2n-row CHP layout we keep JUST the n stabilizer
    generators (x[i], z[i] bitmasks, r[i] sign) and replace the
    destabilizer bookkeeping with an explicit Gaussian-elimination
    membership test for deterministic measurement: Z_w is measured
    deterministically iff it lies in the group spanned by the rows, and
    the outcome is the sign of that expression.  Fewer invariants to
    maintain - nothing to silently break."""
    _scratch = None

    def __init__(self, n: int):
        self.n = n
        self.x = [0] * n
        self.z = [(1 << i) for i in range(n)]
        self.r = [0] * n

    def copy(self) -> "StabState":
        s = StabState.__new__(StabState)
        s.n = self.n
        s.x = list(self.x)
        s.z = list(self.z)
        s.r = list(self.r)
        return s

    # -- gate rules ---------------------------------------------------------
    def _pauli(self, w: int, kind: str) -> None:
        # X_a anticommutes with rows carrying Z_a (and vice versa):
        # the sign of such rows flips
        for i in range(self.n):
            if kind in ("x", "y"):
                self.r[i] ^= (self.z[i] >> w) & 1
            if kind in ("z", "y"):
                self.r[i] ^= (self.x[i] >> w) & 1

    def _h(self, w: int) -> None:
        m = 1 << w
        for i in range(self.n):
            xi, zi = self.x[i], self.z[i]
            self.r[i] ^= ((xi & zi) >> w) & 1
            self.x[i] = (xi & ~m) | (zi & m)
            self.z[i] = (zi & ~m) | (xi & m)

    def _s(self, w: int) -> None:
        m = 1 << w
        for i in range(self.n):
            self.r[i] ^= ((self.x[i] & self.z[i]) >> w) & 1
            self.z[i] ^= self.x[i] & m

    def _cx(self, c: int, t: int) -> None:
        mc, mt = 1 << c, 1 << t
        for i in range(self.n):
            xc = (self.x[i] >> c) & 1
            zt = (self.z[i] >> t) & 1
            xt = (self.x[i] >> t) & 1
            zc = (self.z[i] >> c) & 1
            self.r[i] ^= xc & zt & (xt ^ zc ^ 1)
            self.x[i] ^= xc << t
            self.z[i] ^= zt << c

    def _cz(self, c: int, t: int) -> None:
        # CZ is symmetric: EACH X picks up the other wire's Z
        # (X_c -> X_c Z_t AND X_t -> X_t Z_c)
        for i in range(self.n):
            self.z[i] ^= ((self.x[i] >> c) & 1) << t
            self.z[i] ^= ((self.x[i] >> t) & 1) << c

    def _swap(self, a: int, b: int) -> None:
        if a > b:
            a, b = b, a      # swap is symmetric; keep the shift positive
        ma, mb = (1 << a), (1 << b)
        m = ma | mb
        shift = b - a
        for i in range(self.n):
            xa, xb = self.x[i] & ma, self.x[i] & mb
            self.x[i] = (self.x[i] & ~m) | (xa << shift) | (xb >> shift)
            za, zb = self.z[i] & ma, self.z[i] & mb
            self.z[i] = (self.z[i] & ~m) | (za << shift) | (zb >> shift)

    def _generic_1q(self, w: int, mat) -> None:
        """Numeric conjugation: U P U^dag = i^e P' for each Pauli type;
        phase e folds into the sign bit (mod 2)."""
        table = _conj_table(mat)
        m = 1 << w
        for i in range(self.n):
            xb = (self.x[i] >> w) & 1
            zb = (self.z[i] >> w) & 1
            e, xp, zp = table[(xb, zb)]
            self.r[i] ^= e & 1
            self.x[i] = (self.x[i] & ~m) | (xp << w)
            self.z[i] = (self.z[i] & ~m) | (zp << w)

    # -- public API ----------------------------------------------------------
    def apply(self, g: Gate) -> None:
        n = g.name
        if n == "cx":
            self._cx(g.qubits[0], g.qubits[1])
        elif n == "cz":
            self._cz(g.qubits[0], g.qubits[1])
        elif n == "swap":
            self._swap(g.qubits[0], g.qubits[1])
        elif n in ("x", "y", "z"):
            self._pauli(g.qubits[0], n)
        elif n == "h":
            self._h(g.qubits[0])
        elif n == "s":
            self._s(g.qubits[0])
        elif n == "sdg":
            self._s(g.qubits[0])
            self._pauli(g.qubits[0], "z")
        elif n == "sx":
            self._h(g.qubits[0])
            self._s(g.qubits[0])
            self._h(g.qubits[0])
        elif n in ("rx", "ry", "rz", "p"):
            # generic numeric conjugation: the table construction itself
            # rejects non-Clifford angles, so nothing can slip through
            from .linalg import gate_matrix
            mat = gate_matrix(g)
            self._generic_1q(g.qubits[0], tuple(complex(x) for x in mat))
        else:
            raise ValueError(f"gate {n} is not Clifford")

    def measure(self, w: int, rng: random.Random) -> int:
        """Measure Z on wire w; collapses the state."""
        n = self.n
        anticommuting = [i for i in range(n) if (self.x[i] >> w) & 1]
        if anticommuting:
            # random outcome: make row p the unique X_w-bearing generator,
            # then replace it with the measured Z_w
            p = anticommuting[0]
            for i in anticommuting[1:]:
                self._rowsum(i, p)
            self.x[p] = 0
            self.z[p] = 1 << w
            self.r[p] = rng.getrandbits(1)
            return self.r[p]
        # deterministic: express Z_w in the group via GF(2) elimination
        basis = {}   # pivot column -> (vector, combination mask)

        def insert(v, combo):
            while v:
                c = (v & -v).bit_length() - 1
                if c not in basis:
                    basis[c] = (v, combo)
                    return
                bv, bc = basis[c]
                v ^= bv
                combo ^= bc

        for i in range(n):
            insert(self.x[i] | (self.z[i] << n), 1 << i)
        t = 1 << (n + w)
        combo = 0
        while t:
            c = (t & -t).bit_length() - 1
            if c not in basis:
                raise RuntimeError("Z_w not in stabilizer group")
            bv, bc = basis[c]
            t ^= bv
            combo ^= bc
        rs, xs, zs = 0, 0, 1 << w
        for i in range(n):
            if (combo >> i) & 1:
                rs, xs, zs = self._rowsum_into(rs, xs, zs, self.r[i],
                                               self.x[i], self.z[i])
        assert xs == 0 and zs == 0, "elimination failed to hit identity"
        return rs

    def _rowsum(self, h: int, i: int) -> None:
        s2 = 2 * (self.r[h] + self.r[i])             + bin(self.z[h] & self.x[i]).count("1")             - bin(self.x[h] & self.z[i]).count("1")
        self.r[h] = (s2 % 4) // 2
        self.x[h] ^= self.x[i]
        self.z[h] ^= self.z[i]

    @staticmethod
    def _rowsum_into(rh, xh, zh, ri, xi, zi):
        s2 = 2 * (rh + ri) + bin(zh & xi).count("1")             - bin(xh & zi).count("1")
        return (s2 % 4) // 2, xh ^ xi, zh ^ zi


def _conj_table(mat):
    """U P U^dag = i^e P' for the four Pauli types, from the 2x2 `mat`
    (flat tuple a,b,c,d).  Returns {(xb, zb): (e, x', z')}."""
    from .linalg import mmul
    a, b, c, d = mat
    U = (a, b, c, d)
    Ud = (a.conjugate(), c.conjugate(), b.conjugate(), d.conjugate())
    paulis = {(0, 0): (1 + 0j, 0j, 0j, 1 + 0j),
              (1, 0): (0j, 1 + 0j, 1 + 0j, 0j),
              (0, 1): (1 + 0j, 0j, 0j, -1 + 0j),
              (1, 1): (0j, -1j, 1j, 0j)}
    want = {(0, 0): (1 + 0j, 0j, 0j, 1 + 0j),
            (1, 0): (0j, 1 + 0j, 1 + 0j, 0j),
            (0, 1): (1 + 0j, 0j, 0j, -1 + 0j),
            (1, 1): (0j, -1j, 1j, 0j)}
    table = {}
    for (xb, zb), P in paulis.items():
        UP = mmul(U, P)
        UP_Ud = mmul(UP, Ud)
        best = None
        for (xp, zp), W in want.items():
            tr = (UP_Ud[0].conjugate() * W[0] + UP_Ud[1].conjugate() * W[1]
                  + UP_Ud[2].conjugate() * W[2] + UP_Ud[3].conjugate() * W[3])
            if abs(abs(tr) / 2 - 1.0) < 1e-9:   # Tr-normalized: |<W,M>| = 1 for Pauli match
                phase = tr / abs(tr)
                e = {(1 + 0j): 0, 1j: 1, (-1 + 0j): 2, (-1j): 3}[
                    complex(round(phase.real), round(phase.imag))]
                best = (e, xp, zp)
                break
        if best is None:
            raise ValueError("gate is not Clifford")
        table[(xb, zb)] = best
    return table


def stab_run_shot(circ: Circuit, rng: random.Random,
                  error_rate: float = 0.0,
                  two_q_error: float | None = None,
                  noise=None) -> str:
    """One noiseless-or-stochastic trajectory; returns the outcome string
    (char i = wire i).  Pauli errors are sampled per gate: each involved
    wire errors independently with probability error_rate (1q) /
    two_q_error (2q, default error_rate).  When a `noise` NoiseModel is
    given, per-gate calibrations (including directional per-pair entries)
    take priority over the scalars."""
    st = StabState(circ.num_qubits)
    for g in circ.ops:
        st.apply(g)
        if noise is not None:
            # per-gate calibration rates (honors per-pair entries)
            rate = noise.gate_error(g)
        else:
            rate = error_rate
            if len(g.qubits) == 2:
                rate = two_q_error if two_q_error is not None else error_rate
        if rate > 0:
            for w in g.qubits:
                if rng.random() < rate:
                    st._pauli(w, rng.choice(("x", "y", "z")))
    bits = []
    for w in range(circ.num_qubits):
        bits.append(str(st.measure(w, rng)))
    return "".join(bits)


def stab_sample(circ: Circuit, shots: int = 1000, seed: int = 0,
                error_rate: float = 0.0,
                two_q_error: float | None = None) -> dict:
    """Sample `shots` trajectories; counts by outcome string."""
    rng = random.Random(seed)
    counts: dict = {}
    for _ in range(shots):
        s = stab_run_shot(circ, rng, error_rate, two_q_error)
        counts[s] = counts.get(s, 0) + 1
    return counts
