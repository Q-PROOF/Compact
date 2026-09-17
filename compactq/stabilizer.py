"""Clifford stabilizer tableaux (CHP formalism) for exact Clifford reasoning.

A Clifford gate U is characterized by its conjugation action on the Pauli
group:  U P U^dag = ±P'.  We track the images of the 2n generator Paulis
X_0..X_{n-1}, Z_0..Z_{n-1} as rows of the Aaronson-Gottesman form

    row = (x, z, e)   representing the operator  i^e . (X^x . Z^z)

with x, z bit-vectors of length n and a phase exponent e mod 4.  This is
operator-exact (rows with x=z=1 hold \mp i Y factors in e), which is what
the equality contract needs; plain stabilizer *simulation* only needs the
sign bit, but synthesis requires the exact images.

Two Clifford circuits are equal up to global phase iff their tableaux are
identical — that is the exactness contract used by `clifford_equal`.

Design note: single-qubit gates are conjugated NUMERICALLY (build the 2x2,
conjugate the four single-qubit Pauli products, match against the 8 signed
Pauli candidates) instead of via hand-derived update rules, so a wrong sign
convention cannot ship silently; `from_circuit` returns None for
non-Clifford input.  CX uses the standard bit-rule (oracle-validated).
"""
from __future__ import annotations

from .circuit import Circuit
from .linalg import gate_matrix, mmul

# 2x2 matrices as flat tuples (a, b, c, d)
_X2 = (0j, 1 + 0j, 1 + 0j, 0j)
_Z2 = (1 + 0j, 0j, 0j, -1 + 0j)

_PAULI_BY_XZ = {(0, 0): (1 + 0j, 0j, 0j, 1 + 0j),
                (1, 0): _X2,
                (0, 1): _Z2,
                (1, 1): mmul(_X2, _Z2)}


def _scale(mat, c):
    return tuple(x * c for x in mat)


def _conj(g_mat, op):
    return mmul(mmul(g_mat, op), _dag(g_mat))


def _dag(m):
    a, b, c, d = m
    return (a.conjugate(), c.conjugate(), b.conjugate(), d.conjugate())


def _pauli_conj_table(g_mat, tol=1e-8):
    """Table {(x, z): (r', x', z')} with G (X^x Z^z) G^dag = (-1)^r' X^x' Z^z'.

    Returns None when G is not Clifford (some image is not a signed Pauli).
    """
    table = {}
    for x in (0, 1):
        for z in (0, 1):
            op = _PAULI_BY_XZ[(x, z)]
            img = _conj(g_mat, op)
            found = None
            for e2 in range(4):
                for x2 in (0, 1):
                    for z2 in (0, 1):
                        cand = _scale(_PAULI_BY_XZ[(x2, z2)], 1j ** e2)
                        dev = max(abs(img[k] - cand[k]) for k in range(4))
                        if dev < tol:
                            found = (e2, x2, z2)
                            break
                    if found:
                        break
                if found:
                    break
            if found is None:
                return None
            table[(x, z)] = found
    return table


class Tableau:
    """Stabilizer tableau: images of the 2n Pauli generators."""

    def __init__(self, num_qubits: int):
        self.n = num_qubits
        # rows 0..n-1: images of X_q; rows n..2n-1: images of Z_q
        self.x = [[0] * self.n for _ in range(2 * self.n)]
        self.z = [[0] * self.n for _ in range(2 * self.n)]
        self.e = [0] * (2 * self.n)
        for q in range(self.n):
            self.x[q][q] = 1            # X_q -> X_q
            self.z[self.n + q][q] = 1   # Z_q -> Z_q

    # ------------------------------------------------------------------
    def apply_cx(self, c: int, t: int):
        """Conjugate every row by CX(c, t)."""
        x, z, e = self.x, self.z, self.e
        for i in range(2 * self.n):
            # In the mod-4 product convention the bit updates alone are the
            # exact conjugation (verified exhaustively against matrix truth
            # over all 64 phase states, both CX directions - no sign term;
            # the familiar AG sign-flip is an artifact of mod-2 phases).
            x[i][t] ^= x[i][c]
            z[i][c] ^= z[i][t]

    def apply_1q_table(self, q: int, table: dict):
        """Conjugate every row by a 1q gate with a precomputed Pauli table."""
        x, z, e = self.x, self.z, self.e
        for i in range(2 * self.n):
            e2, x2, z2 = table[(x[i][q], z[i][q])]
            e[i] = (e[i] + e2) & 3
            x[i][q] = x2
            z[i][q] = z2

    # ------------------------------------------------------------------
    def apply_swap(self, a: int, b: int):
        """Conjugate every row by SWAP(a, b): exchange the wires' columns."""
        x, z = self.x, self.z
        for i in range(2 * self.n):
            x[i][a], x[i][b] = x[i][b], x[i][a]
            z[i][a], z[i][b] = z[i][b], z[i][a]

    def copy(self) -> "Tableau":
        out = Tableau.__new__(Tableau)
        out.n = self.n
        out.x = [row[:] for row in self.x]
        out.z = [row[:] for row in self.z]
        out.e = self.e[:]
        return out

    def key(self):
        """Hashable canonical form (equality up to global phase of U)."""
        return (tuple(tuple(r) for r in self.x),
                tuple(tuple(r) for r in self.z),
                tuple(self.e))


def _gate_1q_table(g, cache={}):
    """(cached) Pauli conjugation table for a 1q gate, or None if not Clifford."""
    key = (g.name, g.params)
    if key in cache:
        return cache[key]
    try:
        mat = gate_matrix(g)
    except Exception:
        cache[key] = None
        return None
    table = _pauli_conj_table(mat)
    cache[key] = table
    return table


def tableau_from_circuit(circ: Circuit):
    """Tableau of a Clifford circuit, or None if any gate is non-Clifford."""
    tab = Tableau(circ.num_qubits)
    for g in circ.ops:
        if len(g.qubits) == 1:
            table = _gate_1q_table(g)
            if table is None:
                return None
            tab.apply_1q_table(g.qubits[0], table)
        elif g.name == "cx" and len(g.qubits) == 2:
            tab.apply_cx(g.qubits[0], g.qubits[1])
        else:
            return None  # non-Clifford 2q gate (cz/swap are handled below)
    return tab


def is_clifford(circ: Circuit) -> bool:
    return tableau_from_circuit(circ) is not None


def clifford_equal(circ_a: Circuit, circ_b: Circuit) -> bool:
    """True iff both are Clifford and equal up to global phase."""
    if circ_a.num_qubits != circ_b.num_qubits:
        raise ValueError("qubit count mismatch")
    ta = tableau_from_circuit(circ_a)
    tb = tableau_from_circuit(circ_b)
    if ta is None or tb is None:
        raise ValueError("clifford_equal requires Clifford circuits")
    return ta.key() == tb.key()


def _expand_cz_swap(circ: Circuit) -> Circuit:
    """Rewrite cz/swap into cx (the tableau engine's native 2q gate)."""
    from .circuit import Gate
    out = []
    for g in circ.ops:
        if g.name == "cz":
            a, b = g.qubits
            out.append(Gate("h", (), (b,)))
            out.append(Gate("cx", (), (a, b)))
            out.append(Gate("h", (), (b,)))
        elif g.name == "swap":
            a, b = g.qubits
            out.append(Gate("cx", (), (a, b)))
            out.append(Gate("cx", (), (b, a)))
            out.append(Gate("cx", (), (a, b)))
        else:
            out.append(g)
    return Circuit(circ.num_qubits, out)


def _expand_cz_swap(circ: Circuit) -> Circuit:
    """Rewrite cz/swap into cx (the tableau engine's native 2q gate)."""
    from .circuit import Gate
    out = []
    for g in circ.ops:
        if g.name == "cz":
            a, b = g.qubits
            out.append(Gate("h", (), (b,)))
            out.append(Gate("cx", (), (a, b)))
            out.append(Gate("h", (), (b,)))
        elif g.name == "swap":
            a, b = g.qubits
            out.append(Gate("cx", (), (a, b)))
            out.append(Gate("cx", (), (b, a)))
            out.append(Gate("cx", (), (a, b)))
        else:
            out.append(g)
    return Circuit(circ.num_qubits, out)


def tableau(circ: Circuit):
    """Tableau of a Clifford circuit (cz/swap supported), else None."""
    return tableau_from_circuit(_expand_cz_swap(circ))


def is_clifford(circ: Circuit) -> bool:
    """True iff every gate of the circuit is a Clifford gate."""
    return tableau(circ) is not None


def clifford_equal(circ_a: Circuit, circ_b: Circuit) -> bool:
    """True iff both circuits are Clifford and equal up to global phase.

    Raises ValueError when either circuit is not Clifford or the qubit
    counts differ — callers wanting a plain boolean should gate on
    `is_clifford` first.
    """
    if circ_a.num_qubits != circ_b.num_qubits:
        raise ValueError("qubit count mismatch")
    ta = tableau(circ_a)
    tb = tableau(circ_b)
    if ta is None or tb is None:
        raise ValueError("clifford_equal requires Clifford circuits")
    return ta.key() == tb.key()
