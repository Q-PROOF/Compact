# compactq-check

**Independent checker for Compact compilation certificates.**

Given an original circuit, an optimized circuit, and the certificate that
Compact emitted, this program re-derives the equivalence claim from
scratch — rebuilding unitaries, Clifford tableaux, or phase-polynomial
tables per the certificate's witness kind — and prints a verdict.

Design contract:

- **Zero dependencies.** Pure Python ≥ 3.9.
- **Shares no code with compactq.** Enforced by tests.
- **Never trusts the certificate**: every hash is recomputed, every
  circuit is re-parsed, every equivalence is re-proven from gate lists.
- Small enough to audit in an afternoon.

## Usage

```bash
compactq-check original.qasm optimized.qasm certificate.json
# exit 0 = VALID, 1 = INVALID, 2 = MALFORMED
compactq-check original.qasm optimized.qasm certificate.json --json
```
