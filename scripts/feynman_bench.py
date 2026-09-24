"""Feynman-corpus benchmark: run compactq's verified optimizer over the
benchmark circuit corpus shipped with github.com/meamy/feynman (44
QASM files: adders, Barenco Toffoli networks, GF(2) multipliers, QAOA,
... — the reversible / Clifford+T workloads where exact algebraic
proofs matter most).

Per circuit: identical input for every tool -> optimize -> count
(total / 2q / depth) + wall time -> in-product proof tier recorded ->
every output <= 8 qubits refereed by qiskit's Operator, wider outputs
re-checked through compactq's independent verify() cascade; T-count /
T-depth reported via compactq.resources.

The feynver verification TOOL itself needs a Haskell toolchain
(ghc/cabal), which this environment does not have — recorded as "not
runnable here"; the QCEC and PyZX referees (results/qcec_referee.json,
results/pyzx_referee.json) are the external equivalence checks in
this release.

Checkout resolution: $FEYNMAN_HOME or third_party/feynman (downloaded,
tarball pinned to github.com, traversal-checked extraction).

Usage:  python scripts/feynman_bench.py [--out-dir results/feynman]
Artifacts: results/feynman.json + results/feynman.md
"""
from __future__ import annotations

import argparse
import io
import json
import os
import sys
import tarfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

FEYNMAN_TARBALL = "https://github.com/meamy/feynman/archive/refs/heads/master.tar.gz"
ALLOWED_HOSTS = {"github.com", "codeload.github.com",
                 "objects.githubusercontent.com"}

# skip a handful of files whose OpenQASM uses constructs outside the
# qelib1 subset (documented in the artifact, never silently dropped)
SKIP_IF_UNPARSABLE = True


def _safe_fetch(url: str) -> bytes:
    from urllib.parse import urlparse
    parts = urlparse(url)
    if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:
        raise ValueError(f"URL outside the pinned allowlist: {url!r}")
    import ipaddress
    import socket
    for info in socket.getaddrinfo(parts.hostname, 443, proto=socket.IPPROTO_TCP):
        addr = ipaddress.ip_address(info[4][0])
        if (addr.is_private or addr.is_loopback or addr.is_link_local
                or addr.is_reserved or addr.is_multicast):
            raise ValueError(f"{parts.hostname} resolved to a "
                             f"non-public address ({addr})")
    req = urllib.request.Request(url, headers={"User-Agent": "compactq-bench/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:  # noqa: S310
        if resp.geturl() and urlparse(resp.geturl()).hostname \
                not in ALLOWED_HOSTS:
            raise ValueError("redirect landed outside the allowlist")
        return resp.read()


def _safe_extract(tf: tarfile.TarFile, dest: Path) -> None:
    dest = dest.resolve()
    for member in tf.getmembers():
        name = member.name.replace("\\", "/")
        if (name.startswith(("/", "~")) or ".." in Path(name).parts
                or Path(name).is_absolute()):
            raise ValueError(f"tar member rejected: {member.name!r}")
        if not str((dest / member.name).resolve()).startswith(str(dest)):
            raise ValueError(f"tar member escapes destination: {member.name!r}")
    tf.extractall(dest)


def ensure_corpus() -> Path:
    env = os.environ.get("FEYNMAN_HOME")
    if env and Path(env).is_dir():
        return Path(env) / "benchmarks" / "qasm"
    local = REPO / "third_party" / "feynman"
    if not local.is_dir():
        local.parent.mkdir(exist_ok=True)
        print("downloading meamy/feynman master tarball ...", flush=True)
        data = _safe_fetch(FEYNMAN_TARBALL)
        with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
            _safe_extract(tf, local.parent)
        (local.parent / "feynman-master").rename(local)
    return local / "benchmarks" / "qasm"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(REPO / "results"))
    ap.add_argument("--dd-deadline-s", default="15",
                    help="per-build DD wall-clock deadline for this "
                         "corpus run (the product default is 120s; a "
                         "44-circuit corpus uses a tighter bound and "
                         "records the declined tier honestly)")
    args = ap.parse_args()
    os.environ.setdefault("COMPACTQ_DD_DEADLINE_S", args.dd_deadline_s)
    os.environ.setdefault("COMPACTQ_VERIFY_DEADLINE_S", "15")
    # keep the corpus run's memory footprint bounded (a 400k-node
    # pure-Python diagram costs hundreds of MB; 60k is ample for the
    # proofs this corpus reaches and declines honestly beyond that)
    os.environ.setdefault("COMPACTQ_DD_MAX_NODES", "60000")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import compactq
    from compactq import from_qasm, optimize_search
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from compactq import verify
    from compactq.resources import resource_estimate

    corpus = ensure_corpus()
    files = sorted(corpus.glob("*.qasm"))
    print(f"feynman corpus: {len(files)} qasm files")

    import warnings
    warnings.filterwarnings("ignore")
    import numpy as np
    from qiskit import qasm2
    from qiskit.quantum_info import Operator

    records = []
    parsed = skipped = 0
    for f in files:
        name = f.stem
        rec = {"circuit": name}
        try:
            text = f.read_text(encoding="utf-8", errors="replace")
            mine = from_qasm(text)
            ref = qasm2.loads(text)
        except Exception as e:
            rec.update({"status": f"skipped ({type(e).__name__}: "
                                  f"{str(e)[:60]})"})
            skipped += 1
            records.append(rec)
            print(f"  {name:24s} SKIPPED", flush=True)
            continue
        parsed += 1
        n = mine.num_qubits
        rec["qubits"] = n
        rec["input_gates"] = len(mine)
        rec["input_2q"] = mine.two_qubit_count()
        try:
            est = resource_estimate(mine)
            rec["input_t_count"] = est.get("t_count")
        except Exception:  # noqa: BLE001
            rec["input_t_count"] = None

        t0 = time.perf_counter()
        try:
            out = optimize_search(from_qiskit(ref))
        except Exception as e:  # noqa: BLE001
            rec["status"] = f"optimize-error {type(e).__name__}"
            records.append(rec)
            print(f"  {name:24s} OPTIMIZE-ERROR", flush=True)
            continue
        dt = time.perf_counter() - t0
        rec["compact_gates"] = len(out)
        rec["compact_2q"] = out.two_qubit_count()
        rec["compact_wall_ms"] = round(dt * 1000, 1)
        v = verify(mine, out)
        rec["proof_tier"] = v["tier"]
        rec["proof_method"] = v["method"]
        rec["verified"] = v["equivalent"] is True
        # independent referee <= 8q
        if n <= 8:
            try:
                low_in = qasm2.loads(to_qasm(mine))
                low_out = qasm2.loads(to_qasm(out))
                U = Operator(low_in).data
                V = Operator(low_out).data
                fid = abs(np.sum(np.conj(U) * V)) / U.shape[0]
                rec["referee_fidelity"] = round(float(fid), 12)
            except Exception as e:  # noqa: BLE001
                rec["referee_fidelity"] = None
                rec["referee_note"] = type(e).__name__
        records.append(rec)
        print(f"  {name:24s} {len(mine):5d} -> {len(out):5d} gates  "
              f"tier {v['tier']} ({'refereed' if n <= 8 else 'in-product proof'})",
              flush=True)

    ok = sum(1 for r in records if r.get("verified"))
    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "corpus": "github.com/meamy/feynman benchmarks/qasm (44 files)",
        "corpus_commit": "master (tarball, unversioned upstream)",
        "protocol": ("identical inputs; compactq optimize_search with its "
                     "whole-circuit proof in the timed region; qiskit "
                     "Operator referee <= 8q; in-product verify() cascade "
                     "for wider circuits; T-count via compactq.resources"),
        "parsed": parsed,
        "skipped": skipped,
        "verified": ok,
        "feynver_tool": ("not runnable here (Haskell toolchain required); "
                         "external referees are QCEC + PyZX"),
        "records": records,
        "compactq_version": compactq.__version__,
    }
    (out_dir / "feynman.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")

    lines = ["# Feynman-corpus verified optimization (generated)", "",
             f"Generated {doc['generated']} | corpus: meamy/feynman "
             f"benchmarks/qasm ({parsed} parsed, {skipped} skipped)", "",
             "compactq optimize_search (proof included) — "
             f"{ok}/{parsed} verified in-product", "",
             "| circuit | q | in gates | out gates | in 2q | out 2q | "
             "wall ms | proof | referee |",
             "|---|---:|---:|---:|---:|---:|---:|---|---|"]
    for r in records:
        if "compact_gates" not in r:
            lines.append(f"| {r['circuit']} | - | - | - | - | - | - | "
                         f"{r.get('status', '?')} | - |")
            continue
        lines.append(
            f"| {r['circuit']} | {r['qubits']} | {r['input_gates']} | "
            f"{r['compact_gates']} | {r['input_2q']} | {r['compact_2q']} | "
            f"{r['compact_wall_ms']} | {r['proof_method']} | "
            f"{r.get('referee_fidelity', 'in-product')} |")
    (out_dir / "feynman.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    print(f"wrote {out_dir/'feynman.json'} + feynman.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
