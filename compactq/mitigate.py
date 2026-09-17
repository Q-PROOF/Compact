"""Measurement error mitigation (readout calibration + inversion).

Per-wire readout confusion is tensored: the full 2^n confusion is the
Kronecker product of per-wire 2x2 matrices, so its inverse applies one
wire at a time in O(n * 2^n).  Output is a quasi-distribution: values may
be slightly negative (that is what honest inversion looks like); callers
can clip+renormalize for probabilities.

- calibration_circuits(n): one X-prep circuit per wire, to measure the
  confusion on hardware (feed the resulting counts to calibrate).
- mitigate_counts(counts, calib): corrected quasi-probabilities.
"""
from __future__ import annotations

from .circuit import Circuit, Gate
from .noise import NoiseModel, default_model


def calibration_circuits(n: int) -> list[Circuit]:
    """n circuits: identity (all-0 reference) + X on each wire (all-1)."""
    circs = [Circuit(n, [])]
    for w in range(n):
        circs.append(Circuit(n, [Gate("x", (), (w,))]))
    return circs


def calibrate_from_counts(n: int, counts_by_circuit: list[dict]) -> dict:
    """Derive per-wire (P(1|0), P(0|1)) from calibration-circuit counts.

    counts_by_circuit: [all0_counts, x_on_0_counts, x_on_1_counts, ...],
    each mapping bitstring -> count.  Bitstring convention: character i is
    wire i (little-endian like the rest of compactq).
    """
    assert len(counts_by_circuit) == n + 1
    readout = {}
    base = counts_by_circuit[0]
    tot0 = sum(base.values()) or 1
    for w in range(n):
        xc = counts_by_circuit[1 + w]
        tot1 = sum(xc.values()) or 1
        # P(1|0): fraction of 1s on wire w in the all-0 prep
        p10 = sum(v for k, v in base.items()
                  if k[w] == "1") / tot0
        # P(0|1): fraction of 0s on wire w in the X-on-w prep
        p01 = sum(v for k, v in xc.items()
                  if k[w] == "0") / tot1
        readout[w] = (min(max(p10, 0.0), 0.5), min(max(p01, 0.0), 0.5))
    return readout


def _inverse_2x2(p10, p01):
    """Inverse of [[1-p10, p01], [p10, 1-p01]] (rows = true, cols = observed)."""
    det = (1 - p10) * (1 - p01) - p10 * p01
    if abs(det) < 1e-12:
        return ((1.0, 0.0), (0.0, 1.0))
    return (((1 - p01) / det, -p01 / det),
            (-p10 / det, (1 - p10) / det))


def mitigate_counts(counts: dict, readout: dict, n: int,
                    clip: bool = True) -> dict:
    """Apply the tensored readout inverse to a counts distribution.

    counts: bitstring (char i = wire i) -> count.  Returns bitstring ->
    quasi-probability (clipped+renormalized when clip=True)."""
    tot = sum(counts.values()) or 1
    probs = {k: v / tot for k, v in counts.items()}
    for w in range(n):
        inv = _inverse_2x2(*readout.get(w, (0.0, 0.0)))
        new = {}
        for s, p in probs.items():
            bit = s[w]
            row = 1 if bit == "1" else 0
            # observed -> true contributions on both values of wire w
            for nb, coef in ((0, inv[row][0]), (1, inv[row][1])):
                if abs(coef) < 1e-15:
                    continue
                s2 = s[:w] + str(nb) + s[w + 1:]
                new[s2] = new.get(s2, 0.0) + coef * p
        probs = new
    if clip:
        neg = sum(-v for v in probs.values() if v < 0)
        probs = {k: max(v, 0.0) for k, v in probs.items()}
        s = sum(probs.values())
        if s > 0:
            probs = {k: v / s for k, v in probs.items()}
    return probs


def mitigate_with_model(counts: dict, noise: NoiseModel, clip: bool = True) -> dict:
    """Convenience: mitigate using a NoiseModel's readout calibration."""
    n = max((len(k) for k in counts), default=0)
    readout = {w: noise.readout_error(w) for w in range(n)}
    return mitigate_counts(counts, readout, n, clip=clip)


def mitigate_mle(counts: dict, readout: dict, n: int,
                 iters: int = 60, tol: float = 1e-10) -> dict:
    """Maximum-likelihood correction under the tensored readout model.

    Richardson-Lucy / EM iteration: p <- p * A^T (obs / (A p)), where A
    is the tensored per-wire confusion matrix (rows = true, cols =
    observed).  EM provably increases the likelihood every step and keeps
    the estimate non-negative, converging to the constrained MLE - so,
    unlike the matrix inverse, the result is ALWAYS a physical
    distribution (non-negative, summing to 1).

    A is never materialized: both the forward and transpose applications
    fold the per-wire 2x2 blocks in O(n * 2^n).  Pure Python, practical
    to n <= 12.
    """
    size = 1 << n
    tot = sum(counts.values()) or 1
    obs = [0.0] * size
    for s, c in counts.items():
        obs[int(s[::-1][:n], 2)] = c / tot
    uni = 1.0 / size
    p = [0.5 * o + 0.5 * uni for o in obs]   # full-support start

    blocks = [(w,) + readout.get(w, (0.0, 0.0)) for w in range(n)]

    def fold(vec, transpose):
        """Apply the tensored A (or A^T) to a 2^n vector."""
        out = list(vec)
        for (w, p10, p01) in blocks:
            step = 1 << w
            nxt = list(out)
            for base in range(0, size, step * 2):
                for j in range(base, base + step):
                    a, b = out[j], out[j + step]
                    if transpose:
                        # A^T rows are observed: [1-p10, p10; p01, 1-p01]
                        nxt[j] = (1.0 - p10) * a + p10 * b
                        nxt[j + step] = p01 * a + (1.0 - p01) * b
                    else:
                        # A rows are true: [1-p10, p01; p10, 1-p01]
                        nxt[j] = (1.0 - p10) * a + p01 * b
                        nxt[j + step] = p10 * a + (1.0 - p01) * b
            out = nxt
        return out

    for _ in range(iters):
        q = fold(p, transpose=False)
        ratio = [(obs[i] / q[i]) if q[i] > 1e-15 else 0.0
                 for i in range(size)]
        corr = fold(ratio, transpose=True)
        shift = 0.0
        for j in range(size):
            new = p[j] * corr[j]
            shift += abs(new - p[j])
            p[j] = new
        s = sum(p)
        if s > 0:
            p = [x / s for x in p]
        if shift / size < tol:
            break
    return {format(j, "0%db" % n)[::-1]: p[j]
            for j in range(size) if p[j] > 1e-15}
