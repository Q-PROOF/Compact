"""Zero-noise extrapolation (ZNE) - the exact, provable mitigation layer.

Unitary folding scales the noise: the circuit U is replaced by
U (Udag U)^m (global folding, noise scale factor lambda = 2m+1) - an
IDENTITY rewrite (the folded circuit is unitarily identical to the
input, provable by the proof net) whose gate count - and therefore its
stochastic error rate - scales linearly with lambda.  Measuring an
observable at several lambdas and extrapolating to lambda = 0 recovers
a zero-noise estimate WITHOUT any noise model.

Composable with the suppression stack: run the folded circuits through
`suppress_execute`'s machinery (or plain `run_fn`) - suppression
reduces the per-lambda error, ZNE then removes the residual trend.
"""
from __future__ import annotations

from .circuit import Circuit, Gate


def fold_global(circ: Circuit, factor: int) -> Circuit:
    """Global unitary folding at odd scale factor `factor` (1, 3, 5, ...):
    U (Udag U)^m with m = (factor-1)//2.  The result is unitarily
    identical to `circ` - exactly, not approximately."""
    if factor < 1 or factor % 2 == 0:
        raise ValueError("fold factor must be an odd positive integer")
    m = (factor - 1) // 2
    out = list(circ.ops)
    inv = [_dag_gate(g) for g in reversed(circ.ops)]
    for _ in range(m):
        out.extend(inv)
        out.extend(circ.ops)
    return Circuit(circ.num_qubits, out)


def _dag_gate(g: Gate) -> Gate:
    """The inverse of a supported gate ( Clifford gates and rz/p are
    self-inverse or parameter-negating)."""
    if g.name in ("h", "x", "y", "z", "cx", "cz", "swap"):
        return g
    if g.name == "s":
        return Gate("sdg", (), g.qubits)
    if g.name == "sdg":
        return Gate("s", (), g.qubits)
    if g.name in ("rz", "p", "rx", "ry"):
        return Gate(g.name, (-g.params[0],), g.qubits)
    if g.name == "t":
        return Gate("tdg", (), g.qubits)
    if g.name == "tdg":
        return Gate("t", (), g.qubits)
    raise ValueError(f"no inverse rule for gate {g.name}")


def parity_expectation(counts: dict, obs_wires, n: int) -> float:
    """Expectation of the Z-parity observable on `obs_wires`:
    E = sum_s (-1)^{popcount(s & mask)} p_s.  (char i = wire i.)"""
    mask = [0] * n
    for w in obs_wires:
        mask[w] = 1
    tot = sum(counts.values()) or 1
    acc = 0.0
    for s, c in counts.items():
        par = sum(int(s[w]) * mask[w] for w in range(n)) & 1
        acc += (1 if par == 0 else -1) * c
    return acc / tot


def richardson_zero(values, factors):
    """Exact polynomial (Richardson) extrapolation to lambda = 0 through
    all given (lambda, value) points - closed form via Lagrange weights.
    len(values) must be small (typical ZNE: 2-4 points)."""
    pts = [(lam, v) for lam, v in zip(factors, values)]
    total = 0.0
    for i, (li, vi) in enumerate(pts):
        w = 1.0
        for j, (lj, _) in enumerate(pts):
            if i != j:
                w *= lj / (li - lj)
        total += w * vi
    return total


def poly_zero(values, factors, degree: int = 1) -> float:
    """Least-squares polynomial fit in lambda evaluated at 0.  Uses the
    normal equations in O(degree^3 + degree^2 * points) - pure Python."""
    m = degree + 1
    # normal matrix A = V^T V, b = V^T y, with V[i][k] = lam_i^k
    A = [[sum(f ** (k1 + k2) for f in factors) for k2 in range(m)]
         for k1 in range(m)]
    b = [sum(v * f ** k for v, f in zip(values, factors)) for k in range(m)]
    # Gaussian elimination
    aug = [row[:] + [b[i]] for i, row in enumerate(A)]
    n = m
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[piv][col]) < 1e-15:
            return richardson_zero(values, factors)  # degenerate: fall back
        aug[col], aug[piv] = aug[piv], aug[col]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                f = aug[r][col] / aug[col][col]
                for cc in range(col, n + 1):
                    aug[r][cc] -= f * aug[col][cc]
    coef = [aug[i][n] / aug[i][i] for i in range(n)]
    return coef[0]  # fitted value at lambda = 0


def zne_expectation(values, factors, method: str = "richardson") -> float:
    """Extrapolate the measured expectations (one per fold factor) to the
    zero-noise limit.  Methods: 'richardson' (exact polynomial through
    all points), 'linear'/'quadratic' (least squares), 'exponential'
    (two-parameter fit A + B exp(-c lam) via a small scan over c)."""
    if len(values) != len(factors):
        raise ValueError("one value per factor required")
    if method == "richardson":
        return richardson_zero(values, factors)
    if method == "linear":
        return poly_zero(values, factors, degree=1)
    if method == "quadratic":
        return poly_zero(values, factors, degree=2)
    if method == "exponential":
        best = None
        n = len(values)
        for ci in range(1, 61):          # c in a log-spaced scan
            c = ci * 0.05
            # linearize: E - A ~ B e^{-c lam}; fit A = mean, then LS on logs
            # simple two-point-per-c closure: fit B, c by grid, A by mean
            ws = [math_exp(-c * f) for f in factors]
            wmean = sum(ws) / n
            vmean = sum(values) / n
            bnum = sum((w - wmean) * (v - vmean) for w, v in zip(ws, values))
            bden = sum((w - wmean) ** 2 for w in ws)
            if bden < 1e-15:
                continue
            b = bnum / bden
            a = vmean - b * wmean
            resid = sum((a + b * w - v) ** 2 for w, v in zip(ws, values))
            if best is None or resid < best[0]:
                best = (resid, a)
        return best[1] if best else richardson_zero(values, factors)
    raise ValueError(f"unknown ZNE method {method}")


def math_exp(x: float) -> float:
    # local import guard keeps module import cost zero on Python startup
    import math
    return math.exp(x)


def zne_execute(circ: Circuit, obs_wires, noise, run_fn=None, *,
                factors=(1, 3, 5), shots: int = 3000, seed: int = 0,
                method: str = "richardson", eps_coh: float = 0.0) -> dict:
    """One-call ZNE: folds `circ` to each scale factor, executes through
    `run_fn` (or the built-in simulator), measures the Z-parity on
    `obs_wires`, and extrapolates to the zero-noise limit.

    Every folded circuit is an identity rewrite - verify with the proof
    net once and the whole family is certified.  Returns the per-factor
    expectations, the zero-noise estimate, and the unmitigated reference.
    """
    from .simulate import simulate_counts
    noise_eff = noise
    per_factor = []
    raw_counts = None
    for k, f in enumerate(factors):
        folded = fold_global(circ, f)
        if run_fn is not None:
            counts = run_fn(folded, seed=seed + 31 * k, shots=shots)
        else:
            counts = simulate_counts(folded, noise_eff, shots=shots,
                                     seed=seed + 31 * k, eps_coh=eps_coh)
        e = parity_expectation(counts, obs_wires, circ.num_qubits)
        per_factor.append(e)
        if f == 1:
            raw_counts = counts
    zero = zne_expectation(per_factor, list(factors), method=method)
    return {"zero_noise": zero,
            "per_factor": dict(zip(factors, per_factor)),
            "unmitigated": per_factor[0],
            "factors": list(factors),
            "method": method}
