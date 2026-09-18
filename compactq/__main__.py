"""compactq command line entry point (implementation: compactq.cli).

Usage:
  compactq INPUT.qasm [-o OUTPUT.qasm] [--stats] [--objective OBJ] [--json]
  compactq verify ORIGINAL.qasm OPTIMIZED.qasm [-o certificate.json]
"""
from compactq.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
