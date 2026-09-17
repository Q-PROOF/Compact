"""Lightweight symbolic angles for structure-only circuit optimization.

A parameterized angle is a *linear* expression over named parameters:
``Sym`` carries ``{name: coefficient}`` pairs plus a constant.  Every
operation the structural passes need (add, subtract, negate, scale by a
float) stays exact; no trig, no sympy.  Plain floats bypass the class
entirely, so numeric circuits keep the fast path.

    from compactq.symbolic import param, structure_optimize, bind

    t = param("theta")
    circ = ...   # Gate params may hold Sym expressions
    tmpl = structure_optimize(circ)      # angle-agnostic rewrites to fixpoint
    out = bind(tmpl, {"theta": 0.7})     # numeric circuit, ready to run
"""
from __future__ import annotations


class Sym:
    """Linear expression: sum(coeff_i * name_i) + const."""

    __slots__ = ("terms", "const")

    def __init__(self, terms, const=0.0):
        cleaned = {}
        for name, coeff in dict(terms).items():
            coeff = float(coeff)
            if abs(coeff) > 1e-15:
                cleaned[name] = coeff
        self.terms = cleaned
        self.const = float(const)

    # -- construction ------------------------------------------------------
    @staticmethod
    def const_expr(v):
        return Sym({}, v)

    # -- arithmetic --------------------------------------------------------
    def _coerce(self, other):
        if isinstance(other, Sym):
            return other
        if isinstance(other, (int, float)):
            return Sym({}, float(other))
        return None

    def __add__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        terms = dict(self.terms)
        for k, c in o.terms.items():
            terms[k] = terms.get(k, 0.0) + c
        return Sym(terms, self.const + o.const)

    __radd__ = __add__

    def __neg__(self):
        return Sym({k: -c for k, c in self.terms.items()}, -self.const)

    def __sub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self + (-o)

    def __rsub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return o + (-self)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Sym({k: c * other for k, c in self.terms.items()},
                       self.const * other)
        return NotImplemented

    __rmul__ = __mul__

    # -- queries -----------------------------------------------------------
    @property
    def is_const(self):
        return not self.terms

    def value(self, mapping=None):
        """Numeric value; requires all parameters bound when mapping=None."""
        v = self.const
        for k, c in self.terms.items():
            if mapping is None:
                if not isinstance(k, float):
                    raise ValueError(f"unbound parameter {k!r}")
                v += c * k
            else:
                v += c * float(mapping[k])
        return v

    # -- identity ----------------------------------------------------------
    def __eq__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        if self.const != o.const or set(self.terms) != set(o.terms):
            return False
        return all(abs(self.terms[k] - o.terms[k]) < 1e-12 for k in self.terms)

    def __hash__(self):
        return hash((tuple(sorted(self.terms.items())), self.const))

    def __repr__(self):
        parts = [f"{c}*{k}" for k, c in self.terms.items()]
        if self.const or not parts:
            parts.append(repr(self.const))
        return "Sym(" + " + ".join(parts) + ")"


def param(name):
    """A symbolic angle leaf."""
    return Sym({name: 1.0}, 0.0)


def is_zero(x):
    """Exact-ish zero test that works for floats and Sym."""
    if isinstance(x, Sym):
        return x.is_const and abs(x.const) < 1e-12
    return abs(x) < 1e-12


def wrap(x):
    """Wrap numeric angles into (-pi, pi]; symbolic angles pass through
    unwrapped (RZ remains exact for unwrapped angles)."""
    import math
    if isinstance(x, Sym):
        return x
    return (x + math.pi) % (2.0 * math.pi) - math.pi


def bind_gate(g: Gate_like, mapping):
    """Return g with every Sym parameter bound to floats."""
    if not any(isinstance(p, Sym) for p in g.params):
        return g
    new_params = tuple(p.value(mapping) if isinstance(p, Sym) else p
                       for p in g.params)
    return Gate(g.name, new_params, g.qubits)


class _GateLike:  # typing alias only
    pass


Gate_like = object

from .circuit import Gate  # noqa: E402  (after Sym to avoid cycles)


def bind(circ, mapping):
    """Bind every symbolic parameter in `circ` to a float; returns a new
    fully-numeric circuit (the input is untouched)."""
    from .circuit import Circuit
    mapping = {k: float(v) for k, v in mapping.items()}
    return Circuit(circ.num_qubits, [bind_gate(g, mapping) for g in circ.ops])


def structure_optimize(circ, max_iters: int = 12):
    """Angle-agnostic optimization to fixpoint: commutation-driven CX/CZ
    cancellation, the CX-phase rewrite (folds symbolic angles exactly),
    template reorderings and 1q slides.  Excludes every pass that needs
    numeric angles (peephole folding, KAK, parity, cliffordize)."""
    from .transforms import commute_cancel, slide_1q, swap_template
    from .cp_pass import cancel_cx_through_diagonal, merge_cp
    from .templates import template_pass

    def step(c):
        c = commute_cancel(c)
        c = merge_cp(cancel_cx_through_diagonal(c))
        c = template_pass(c)
        c = swap_template(slide_1q(c))
        return c

    cur = circ
    best = (cur.two_qubit_count(), len(cur.ops))
    for _ in range(max_iters):
        nxt = step(cur)
        score = (nxt.two_qubit_count(), len(nxt.ops))
        if score >= best:
            break
        cur, best = nxt, score
    return cur
