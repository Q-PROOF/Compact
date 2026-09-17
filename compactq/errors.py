"""Public exception types for compactq."""


class UnsupportedCircuitError(ValueError):
    """Raised when an input circuit contains operations the unitary
    optimizer cannot represent without changing program semantics:
    mid-circuit measurement, reset, or classical control flow.

    Trailing measurements at the very end of an OpenQASM program are
    dropped (the unitary core is optimized); anything earlier is rejected
    because the optimizer would otherwise silently change what the
    program computes.  The Qiskit bridge rejects measurement and reset
    outright - strip them first, e.g.
    ``qc.remove_final_measurements(inplace=True)``.
    """
