# Security Policy

## Supported versions

| version | support |
|---|---|
| 0.1.x | security fixes |
| < 0.1 | n/a — not publicly released |

## Reporting

Use GitHub **private vulnerability reporting** (Security tab → Report a
vulnerability) rather than a public issue. Include the affected version,
a minimal reproducer, and the proof status the tool reported.

## Scope notes

compactq processes circuit *text and structures* locally; it has no
network code, no telemetry, and no dynamic execution of input files.
The realistic security surface is therefore:

- parser robustness against malformed OpenQASM input (denial of service
  via pathological input, or crashes)
- dependency-free guarantee violations (an unexpected import would be a
  supply-chain regression — CI enforces the zero-dep core)

Reports on either are welcome.
