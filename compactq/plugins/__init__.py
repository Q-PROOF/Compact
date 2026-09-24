"""Optional ecosystem integrations for compactq.

`qiskit_plugin` — Compact as a Qiskit ``TransformationPass`` with the
verification verdict preserved as pass output (extras:
``compactq[qiskit]``).

`pytket_plugin` — the same contract for pytket, usable directly or via
``as_tket_pass()`` inside ``SequencePass``/``CompilationUnit`` workflows
(extras: ``compactq[pytket]``).

Both modules import their toolchain lazily: the compactq core stays
zero-dependency.
"""
