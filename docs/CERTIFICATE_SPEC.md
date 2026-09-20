# Certificate Specification v1 — `whole_circuit_v1`

**Status:** implemented (`compactq.cert`, checked by the independent
`compactq-check` program)
**Composition:** `whole_circuit_v1` — the certificate covers the entire
input→output pair with one witness. Region-level granularity
(`tensor_local_v1`) is specified for the next phase and reserved.

## Guarantees

A certificate asserts: *the circuit serialized in `output.qasm` is
unitarily equivalent to the circuit serialized in `input.qasm`, up to a
global phase, as re-computed by the witness method below.*  It says
nothing about whether the input matches the user's intent, nor about
hardware execution.

If a sound witness cannot be produced for the whole circuit, **no
certificate is issued** — the emitter returns `certificate: None` with a
machine-readable reason.  Certificates are never faked or partially
covered.

## Schema

```jsonc
{
  "cert_version": 1,
  "composition": "whole_circuit_v1",
  "producer": { "name": "compactq", "version": "0.2.2" },
  "input":  { "qasm": "<OpenQASM 2.0>", "sha256": "<hex>" },
  "output": { "qasm": "<OpenQASM 2.0>", "sha256": "<hex>" },
  "n_qubits": 4,
  "global_phase_convention": "equal_up_to_global_phase",
  "witness": {
    "kind": "unitary",          // unitary | stabilizer | phase_polynomial
    "n_qubits": 4,
    "tolerance": 1e-7
  },
  "claims": { "exact": true, "approximate": false, "permuted": false }
}
```

### Witness kinds (strongest applicable is chosen)

| kind | covers | checker re-computation | width limit |
|---|---|---|---|
| `unitary` | any circuit | rebuild both unitaries from the stored gate lists; require \|Tr(A†B)\|/d ≥ 1 − tolerance | ≤ 8 qubits |
| `stabilizer` | Clifford circuits | rebuild Aaronson–Gottesman tableaux for both; require tableau equality (phase-exact up to global phase) | any width |
| `phase_polynomial` | CNOT + diagonal (RZ/P/S/T…) circuits | recompute the parity→angle table for both; require identical parities and angles equal mod 2π within tolerance | any width |

The checker **re-derives every quantity from the stored gate lists** —
witnesses carry no trusted numbers.  Hash mismatches (certificate vs the
two provided QASM files) are a hard INVALID.

## Verdicts

- `VALID` — every re-computation confirms the claim.
- `INVALID` — well-formed certificate, failed equivalence (a soundness
  event: block the release, file a bug).
- `MALFORMED` — unparseable schema/circuits (exit code 2; never a
  soundness signal).
