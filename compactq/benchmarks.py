"""Standard benchmark circuits, generated locally (no downloads needed)."""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate


def ghz(n: int) -> Circuit:
    c = Circuit(n, [])
    c.append(Gate("h", (), (0,)))
    for k in range(1, n):
        c.append(Gate("cx", (), (0, k)))
    return c


def _cp(c: Circuit, theta: float, ctrl: int, tgt: int) -> None:
    """Controlled-phase via the standard 2-CX identity (exact up to phase)."""
    c.append(Gate("p", (theta / 2,), (ctrl,)))
    c.append(Gate("cx", (), (ctrl, tgt)))
    c.append(Gate("p", (-theta / 2,), (tgt,)))
    c.append(Gate("cx", (), (ctrl, tgt)))
    c.append(Gate("p", (theta / 2,), (ctrl,)))


def qft(n: int) -> Circuit:
    c = Circuit(n, [])
    for i in range(n):
        c.append(Gate("h", (), (i,)))
        for j in range(i + 1, n):
            _cp(c, math.pi / (2 ** (j - i)), j, i)
    return c


def clifford_ladder(n: int, seed: int = 7) -> Circuit:
    """Random Clifford circuit — the optimizer should compress these hard."""
    rng = random.Random(seed)
    c = Circuit(n, [])
    for _ in range(3 * n):
        pick = rng.random()
        q = rng.randrange(n)
        if pick < 0.45:
            c.append(Gate(rng.choice(["h", "s", "sdg", "x", "z"]), (), (q,)))
        else:
            a = rng.randrange(n)
            b = (a + 1 + rng.randrange(n - 1)) % n
            c.append(Gate("cx", (), (a, b)))
    return c


def brickwork(n: int, layers: int, seed: int = 11) -> Circuit:
    """Nearest-neighbour 2-qubit brickwork with rotations (hardware-friendly)."""
    rng = random.Random(seed)
    c = Circuit(n, [])
    for _ in range(layers):
        for q in range(n):
            c.append(Gate(rng.choice(["rz", "rx", "ry"]), (rng.uniform(-3, 3),), (q,)))
        for pair in range(0, n - 1, 2):
            c.append(Gate("cx", (), (pair, pair + 1)))
        for q in range(n):
            c.append(Gate("rz", (rng.uniform(-3, 3),), (q,)))
        for pair in range(1, n - 1, 2):
            c.append(Gate("cx", (), (pair, pair + 1)))
    return c


def random_circuit(n: int, depth: int, seed: int = 0) -> Circuit:
    rng = random.Random(seed)
    c = Circuit(n, [])
    for _ in range(depth):
        pick = rng.random()
        if pick < 0.30:
            c.append(Gate(rng.choice(["h", "t", "tdg", "s", "x"]), (), (rng.randrange(n),)))
        elif pick < 0.60:
            c.append(Gate(rng.choice(["rx", "ry", "rz"]), (rng.uniform(-3, 3),), (rng.randrange(n),)))
        elif pick < 0.92:
            a = rng.randrange(n)
            b = rng.randrange(n)
            while b == a:
                b = rng.randrange(n)
            c.append(Gate("cx", (), (a, b)))
        else:
            a = rng.randrange(n)
            b = rng.randrange(n)
            while b == a:
                b = rng.randrange(n)
            c.append(Gate("swap", (), (a, b)))
    return c


def default_suite() -> dict:
    suite = {
        "ghz-5": ghz(5),
        "qft-3": qft(3),
        "qft-4": qft(4),
        "clifford-ladder-4": clifford_ladder(4, seed=7),
        "clifford-ladder-5": clifford_ladder(5, seed=13),
        "brickwork-4x4": brickwork(4, 4, seed=11),
        "random-4q-40": random_circuit(4, 40, seed=1),
        "random-5q-60": random_circuit(5, 60, seed=2),
    }
    return suite


def load_qasm(path: str) -> Circuit:
    """Load a QASM2 file (convenience wrapper)."""
    from .io_qasm import from_qasm
    with open(path) as f:
        return from_qasm(f.read())
