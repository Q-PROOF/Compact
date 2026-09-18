"""`compactq verify ORIGINAL.qasm OPTIMIZED.qasm` — the certificate CLI.

Independently verifies two circuits with `compactq.verify` and emits the
machine-checkable compilation certificate (JSON).  Exit code 0 = equivalent,
1 = not equivalent (or error).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _read_input(path: str) -> str:
    """Read a user-supplied path safely: expand and fully resolve the
    path (no ambiguity), and require an existing regular file."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise SystemExit(f"input file does not exist: {p}")
    return p.read_text(encoding="utf-8")


def _write_out(dest: str, text: str) -> None:
    """Write CLI output safely: expanded, resolved, existing parent."""
    out_path = Path(dest).expanduser().resolve()
    if not out_path.parent.is_dir():
        raise SystemExit(f"--output directory does not exist: "
                         f"{out_path.parent}")
    out_path.write_text(text, encoding="utf-8")
    print(f"compactq: wrote {out_path}", file=sys.stderr)


def verify_main(argv) -> int:
    import json
    from compactq import build_certificate, from_qasm

    ap = argparse.ArgumentParser(
        prog="compactq verify",
        description="independently verify two circuits and emit a "
                    "machine-checkable compilation certificate")
    ap.add_argument("original", help="original OpenQASM 2.0 file")
    ap.add_argument("optimized", help="optimized OpenQASM 2.0 file")
    ap.add_argument("-o", "--output",
                    help="write certificate JSON here (default: stdout)")
    ap.add_argument("--target", default=None,
                    help="optional hardware target name to record")
    ap.add_argument("--estimated-fidelity", type=float, default=None,
                    help="optional estimated hardware fidelity to record")
    args = ap.parse_args(argv)

    ca = from_qasm(_read_input(args.original))
    cb = from_qasm(_read_input(args.optimized))
    cert = build_certificate(ca, cb, target_name=args.target,
                             estimated_fidelity=args.estimated_fidelity)
    text = json.dumps(cert, indent=2) + "\n"
    if args.output:
        _write_out(args.output, text)
    else:
        sys.stdout.write(text)
    eq = cert["verification"]["equivalent"]
    print(f"compactq: verification: {eq} "
          f"(tier {cert['verification']['tier']}, "
          f"{cert['verification']['method']})", file=sys.stderr)
    return 0 if eq is True else 1
