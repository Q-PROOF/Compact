"""compactq command line interface.

Usage:
  compactq INPUT.qasm [-o OUTPUT.qasm] [--stats] [--approx FIDELITY]
             [--objective {2q,depth,gate_count,weighted,latency}]
             [--no-verify] [--native {cz,ecr,iswap}] [--json]
  compactq verify ORIGINAL.qasm OPTIMIZED.qasm [-o cert.json]

Reads an OpenQASM 2.0 file, optimizes it, and writes the result (stdout by
default).  Exit code is 0 on success.

`compactq verify` independently verifies two circuits and emits the
machine-checkable compilation certificate as JSON (exit 0 = equivalent,
1 = not equivalent or error).

Verification policy (reported on stderr and in --json):
  exact-unitary    whole-circuit dense proof (default, within the prover's
                   qubit limit: 8 with the native kernels, 6 without)
  compositional    block-structured wide circuits: per-block dense proofs
                   composed (verification tier 4)
  randomized-exact large circuits: K=32 random-state verification
                   (probabilistically exact; see compactq.verify_large)
  approximate      --approx mode: per-block fidelity guarantees
  unverified       --no-verify, or circuits beyond the randomized prover's
                   reach (> 30 qubits)

--objective selects which metric may never grow: 2q (default; 2-qubit
count, then gates, then depth), depth, gate_count, latency (depth alias),
or weighted (minimize 1.0*2q + 0.1*depth + 0.02*gates).

Trailing measurements in the input are dropped (the unitary core is
optimized); mid-circuit measurement and reset are rejected.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# dispatch mode for circuits beyond the prover's reach: optimize without the
# whole-circuit proof and report the honest "unverified" status
_NO_PROOF = False


def _read_input(path: str) -> str:
    """Read a user-supplied path safely: expand and fully resolve the
    path (no ambiguity), and require an existing regular file."""
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise SystemExit(f"input file does not exist: {p}")
    return p.read_text(encoding="utf-8")


def _write_out(dest: str, text: str) -> None:
    """Write CLI output to the user-specified destination safely: the
    path is expanded and fully resolved (no ambiguity about where the
    bytes land) and its parent directory must already exist."""
    out_path = Path(dest).expanduser().resolve()
    if not out_path.parent.is_dir():
        raise SystemExit(f"--output directory does not exist: "
                         f"{out_path.parent}")
    out_path.write_text(text, encoding="utf-8")
    print(f"compactq: wrote {out_path}", file=sys.stderr)


def _stats_dict(circ):
    return {
        "gates": len(circ.ops),
        "two_qubit": circ.two_qubit_count(),
        "depth": circ.depth(),
    }


def _verify_main(argv):
    """`compactq verify ORIGINAL.qasm OPTIMIZED.qasm [-o cert.json]`:
    independently verify two circuits and emit the compilation
    certificate (JSON).  Exit code 0 = equivalent, 1 = not (or error)."""
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


def _suppress_main(args, circ, to_qasm):
    """--suppress: the automated suppression pipeline, offline."""
    import json
    import compactq
    from compactq import default_model, NoiseModel
    from compactq.suppress import suppress_plan

    if args.noise:
        noise = NoiseModel.from_dict(circ.num_qubits,
                                     json.loads(_read_input(args.noise)))
    else:
        noise = default_model(circ.num_qubits)

    plan = suppress_plan(circ, noise, variants=max(1, args.variants),
                         dd_sequence=args.dd_sequence)
    report = plan["report"]
    out = plan["variants"][0]

    qasm = to_qasm(out)
    if args.output:
        _write_out(args.output, qasm)
    else:
        sys.stdout.write(qasm)

    print(str(report), file=sys.stderr)
    if args.report:
        _write_out(args.report, json.dumps(report.to_dict(), indent=2) + "\n")

    if args.json:
        payload = {
            "version": compactq.__version__,
            "status": "suppressed",
            "num_qubits": circ.num_qubits,
            "variants": plan["num_variants"],
            "dd_sequence": args.dd_sequence,
            "mapping": plan["mapping"],
            "before": _stats_dict(circ),
            "after": _stats_dict(out),
            "report": report.to_dict(),
        }
        sys.stdout.write(json.dumps(payload, indent=2) + "\n")
    return 0


def _optimize_dispatch(circ, approx, verify, objective="2q"):
    """Returns (circuit, proof_status, dense_limit).  Never raises on
    circuit width: wide circuits fall back to unverified optimization
    with an honest status instead of crashing."""
    from compactq.equivalence import _MAX_QUBITS as dense_limit
    if approx is not None:
        from compactq.approximate import approximate
        return approximate(circ, min_fidelity=approx), "approximate", dense_limit
    if not verify:
        from compactq import optimize_search
        return optimize_search(circ, verify=_NO_PROOF,
                               objective=objective), "unverified", dense_limit
    if circ.num_qubits <= dense_limit:
        from compactq import optimize_search
        return optimize_search(circ, verify=True,
                               objective=objective), "exact-unitary", dense_limit
    try:
        from compactq.verify_large import optimize_large
        out, st = optimize_large(circ)
        if st == "exact":
            return out, "randomized-exact", dense_limit
        # randomized verification rejected the optimized circuit: the
        # original is returned unoptimized rather than shipping unproven
        return out, "verification-rejected (original returned)", dense_limit
    except ValueError:
        # beyond the randomized prover's reach: optimize without proof
        from compactq import optimize_search
        return optimize_search(circ, verify=_NO_PROOF,
                               objective=objective), "unverified", dense_limit


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "verify":
        return _verify_main(argv[1:])

    parser = argparse.ArgumentParser(
        prog="compactq",
        description="compactq: the verified quantum circuit optimizer",
    )
    parser.add_argument("input", help="input OpenQASM 2.0 file ('-' = stdin)")
    parser.add_argument("-o", "--output", help="output file (default: stdout)")
    parser.add_argument("--stats", action="store_true",
                        help="print before/after statistics to stderr")
    parser.add_argument("--approx", type=float, default=None, metavar="FIDELITY",
                        help="approximate mode: per-block fidelity floor "
                             "(e.g. 0.99); trades bounded fidelity for fewer "
                             "2-qubit gates")
    parser.add_argument("--objective", default="2q",
                        choices=["2q", "depth", "gate_count", "weighted",
                                 "latency"],
                        help="metric that may never grow: 2q (default), "
                             "depth, gate_count, latency (= depth), or "
                             "weighted (1.0*2q + 0.1*depth + 0.02*gates)")
    parser.add_argument("--no-verify", action="store_true",
                        help="skip the whole-circuit proof (large circuits "
                             "automatically use randomized state verification)")
    parser.add_argument("--native", default=None, choices=["cz", "ecr", "iswap"],
                        help="re-express 2q gates in a machine-native basis")
    parser.add_argument("--json", action="store_true",
                        help="print a machine-readable result summary (JSON) "
                             "to stdout; the QASM circuit is then written "
                             "only when -o is given")
    parser.add_argument("--suppress", action="store_true",
                        help="run the error-suppression pipeline (optimize, "
                             "layout/route, twirl, decouple) and emit the "
                             "suppressed circuit plus a per-stage report")
    parser.add_argument("--noise", default=None, metavar="JSON",
                        help="with --suppress: NoiseModel calibrations "
                             "(T1_us, T2_us, readout, gate_infidelity, "
                             "drift_rate, durations_ns); defaults to a "
                             "representative superconducting device")
    parser.add_argument("--report", default=None, metavar="JSON",
                        help="with --suppress: write the SuppressionReport "
                             "artifact (per-stage proof levels) here")
    parser.add_argument("--dd-sequence", default="auto",
                        choices=["auto", "xy4", "xy8", "xzx", "pdd4"],
                        help="dynamical-decoupling sequence family "
                             "(default: auto)")
    parser.add_argument("--variants", type=int, default=4, metavar="K",
                        help="with --suppress: number of twirled variants "
                             "to build (default 4)")
    args = parser.parse_args(argv)

    import compactq
    from compactq import from_qasm, to_qasm

    text = sys.stdin.read() if args.input == "-" else _read_input(args.input)
    circ = from_qasm(text)

    if args.suppress:
        return _suppress_main(args, circ, to_qasm)

    out, status, dense_limit = _optimize_dispatch(
        circ, args.approx, verify=not args.no_verify,
        objective=args.objective)
    if args.native is not None:
        from compactq.native import rebase
        out = rebase(out, args.native)

    from compactq.verify import proof_status_tier, VERIFICATION_TIERS
    tier = proof_status_tier(status)
    guarantee = VERIFICATION_TIERS[tier]

    result = {
        "version": compactq.__version__,
        "status": status,
        "verification_tier": tier,
        "verification_guarantee": guarantee,
        "num_qubits": circ.num_qubits,
        "dense_proof_limit": dense_limit,
        "native": args.native,
        "objective": args.objective,
        "before": _stats_dict(circ),
        "after": _stats_dict(out),
    }

    if args.json:
        import json
        sys.stdout.write(json.dumps(result, indent=2) + "\n")
        if args.output:
            _write_out(args.output, to_qasm(out))
    else:
        qasm = to_qasm(out)
        if args.output:
            _write_out(args.output, qasm)
        else:
            sys.stdout.write(qasm)

    print(f"compactq: proof status: {status} "
          f"(tier {tier}: {guarantee})", file=sys.stderr)
    if args.stats:
        print(f"before: {circ.stats()}", file=sys.stderr)
        print(f"after : {out.stats()}", file=sys.stderr)
    return 0
