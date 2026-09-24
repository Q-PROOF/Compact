"""pytket pass for compactq (requires `pip install pytket`).

Wraps the compactq optimizer as a pytket pass, usable directly or inside
a pytket ``CompilationUnit`` / standard pass sequence:

    from compactq.plugins.pytket_plugin import CompactTketPass
    tk2 = CompactTketPass().apply(tk)          # direct

    from pytket.passes import SequencePass
    seq = SequencePass([CompactTketPass(), FullPeepholeOptimise()])
    seq.apply(tk)

The proof is preserved as inspectable output, never silently dropped:
``CompactTketPass.last_proof`` carries the ``compactq.verify`` verdict
(tier, method, wall time) of the most recent ``apply``.

Imported lazily — compactq itself stays zero-dependency.
"""
from __future__ import annotations


class CompactTketPass:
    """compactq as a pytket pass.

    Keyword arguments:
      verify=True         run the whole-circuit independent verification
                          after optimizing and record the verdict as
                          proof metadata (the optimization itself is
                          exact either way; False only skips the
                          re-check).
      approximate=False   bounded-fidelity approximate mode instead of
                          the exact search portfolio.
      min_fidelity=0.99   per-block fidelity floor for approximate mode.
      proof_metadata=True keep ``last_proof`` updated on every apply.

    Instances are re-usable; ``apply`` is pure with respect to the
    circuit (returns a new pytket ``Circuit``).
    """

    def __init__(self, verify: bool = True, approximate: bool = False,
                 min_fidelity: float = 0.99, proof_metadata: bool = True):
        self._verify = verify
        self._approx = approximate
        self._min_fid = min_fidelity
        self._proof_metadata = proof_metadata
        self.last_proof = None

    # ------------------------------------------------------------------
    def _to_compactq(self, tk):
        from ..pytket_bridge import from_pytket
        return from_pytket(tk)

    def _to_tket(self, circ):
        from ..pytket_bridge import to_pytket
        return to_pytket(circ)

    def apply(self, tk):
        """Optimize a pytket ``Circuit``; returns a new ``Circuit``."""
        import compactq
        circ = self._to_compactq(tk)
        if self._approx:
            from ..approximate import approximate
            out = approximate(circ, min_fidelity=self._min_fid)
        else:
            from ..search import optimize_search
            out = optimize_search(circ, verify=self._verify)
        out_tk = self._to_tket(out)
        if self._proof_metadata:
            self.last_proof = dict(compactq.verify(circ, out))
        return out_tk

    # pytket BasePass-style convenience so the wrapper can sit inside a
    # CompilationUnit workflow: cu = CompilationUnit(tk); pass.apply(cu)
    def __call__(self, tk):
        return self.apply(tk)

    # ------------------------------------------------------------------
    def as_tket_pass(self, label: str = "compactq"):
        """Return a real pytket ``BasePass`` (a ``CustomPass`` wrapping
        this optimizer) for use inside ``SequencePass`` /
        ``CompilationUnit`` workflows.  The proof verdict is still
        recorded on this instance's ``last_proof`` after each run."""
        from pytket.passes import CustomPass

        def _transform(circ):
            return self.apply(circ)

        return CustomPass(_transform, label=label)
