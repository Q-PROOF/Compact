"""Standard coupling-map topologies for routing / hardware-mapped benchmarks.

Presets (scaled to `n` qubits):
  line        nearest-neighbor chain (n-1 edges)
  grid        square grid, side = ceil(sqrt(n))
  heavy_hex   IBM-style heavy-hex lattice (rows with vertical links on
              even columns; the standard IBM Eagle/Hexagonal pattern,
              scaled to n qubits)
  all_to_all  complete graph (ion-trap style), n(n-1)/2 edges
"""
from __future__ import annotations

import math

__all__ = ["coupling_preset", "PRESETS"]

PRESETS = ("line", "grid", "heavy_hex", "all_to_all")


def coupling_preset(name: str, n: int) -> list[set[int]]:
    """Return the coupling map for a named topology over `n` qubits."""
    n = int(n)
    if n < 2:
        return []
    key = name.lower().replace("-", "_").replace(" ", "_")
    if key == "line":
        return [{j, j + 1} for j in range(n - 1)]
    if key == "all_to_all":
        return [{a, b} for a in range(n) for b in range(a + 1, n)]
    if key == "grid":
        side = max(2, int(math.ceil(math.sqrt(n))))
        pairs = set()
        for r in range(side):
            for c in range(side):
                q = r * side + c
                if q >= n:
                    continue
                if c + 1 < side and q + 1 < n:
                    pairs.add((q, q + 1))
                if r + 1 < side and q + side < n:
                    pairs.add((q, q + side))
        return [set(p) for p in sorted(pairs)]
    if key in ("heavy_hex", "heavyhex"):
        cols = max(3, int(math.ceil(math.sqrt(n))) | 1)  # odd column count
        idx: dict[tuple[int, int], int] = {}
        rows = []
        q = 0
        r = 0
        while q < n:
            row = []
            for c in range(cols):
                if q < n:
                    idx[(r, c)] = q
                    row.append(q)
                    q += 1
            rows.append(row)
            r += 1
        pairs = set()
        for row in rows:
            for i in range(len(row) - 1):
                pairs.add((row[i], row[i + 1]))
        for r in range(len(rows) - 1):
            for (rr, cc), qi in idx.items():
                if rr == r and cc % 2 == 0 and (r + 1, cc) in idx:
                    pairs.add((qi, idx[(r + 1, cc)]))
        return [set(p) for p in sorted(pairs)]
    raise ValueError(f"unknown topology {name!r}; expected one of {PRESETS}")
