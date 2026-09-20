"""`compactq-check` command line: in.qasm out.qasm cert.json

Exit codes: 0 = VALID, 1 = INVALID, 2 = MALFORMED.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _read(path: str) -> str:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        print(f"compactq-check: input file does not exist: {p}",
              file=sys.stderr)
        raise SystemExit(2)
    return p.read_text(encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="compactq-check",
        description="independent checker for Compact compilation "
                    "certificates (shares no code with compactq)")
    ap.add_argument("original", help="original OpenQASM 2.0 file")
    ap.add_argument("optimized", help="optimized OpenQASM 2.0 file")
    ap.add_argument("certificate", help="certificate JSON emitted by "
                                        "compactq")
    ap.add_argument("--json", action="store_true",
                    help="machine output (one JSON verdict object)")
    args = ap.parse_args(argv)

    import compactq_check.checker as checker

    try:
        cert = json.loads(_read(args.certificate))
        verdict = checker.check_certificate(
            cert, _read(args.original), _read(args.optimized))
    except ValueError as e:
        print(json.dumps({"verdict": "MALFORMED", "detail": str(e)}))
        return 2
    except json.JSONDecodeError as e:
        print(json.dumps({"verdict": "MALFORMED",
                          "detail": f"certificate JSON: {e}"}))
        return 2

    if args.json:
        print(json.dumps(verdict, indent=2))
    else:
        line = (f"verdict: {verdict['verdict']} | equivalent: "
                f"{verdict.get('equivalent')} | witness: "
                f"{verdict.get('tier_kind')} ({verdict.get('method')})")
        print(line)
        if verdict["verdict"] != "VALID":
            print(f"detail: {verdict.get('detail', '')}", file=sys.stderr)
    return 0 if verdict["verdict"] == "VALID" else 1


if __name__ == "__main__":
    sys.exit(main())
