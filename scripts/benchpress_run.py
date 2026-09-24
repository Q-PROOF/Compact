"""Run the Qiskit Benchpress abstract-transpilation workout with ONE
gym (compactq or qiskit) per invocation, on the QASMBench suites at the
all-to-all topology, producing pytest-benchmark's native JSON records.

This drives a REAL Benchpress checkout (github.com/qiskit/benchpress):
the compactq gym from `benchpress_integration/` is copied in, the two
upstreamable dispatch elifs are applied, and pytest runs the actual
Benchpress workout module.  Invoke once per gym (a fresh process per
gym keeps Benchpress's singleton Configuration isolated); results land
in results/benchpress/ and are committed unchanged.

Checkout resolution order:
  1. $BENCHPRESS_HOME (an existing clone)
  2. third_party/benchpress (downloaded + extracted on first run;
     gitignored).  Download is pinned to https on github.com only and
     tar members are traversal-checked before extraction.

Usage:
  python scripts/benchpress_run.py --gym compactq --sizes small,medium

Requires: pytest, pytest-benchmark.
"""
from __future__ import annotations

import argparse
import io
import ipaddress
import json
import os
import shutil
import socket
import sys
import tarfile
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
GYM_SRC = REPO / "benchpress_integration" / "compactq_gym"
PATCH = REPO / "benchpress_integration" / "patches" / "apply_patches.py"
BENCHPRESS_TARBALL = "https://github.com/qiskit/benchpress/archive/refs/heads/main.tar.gz"
ALLOWED_HOSTS = {"github.com", "codeload.github.com",
                 "objects.githubusercontent.com", "raw.githubusercontent.com"}
KNOWN_SIZES = {"small", "medium", "large"}


def _safe_fetch(url: str) -> bytes:
    """Fetch with scheme + host allowlist + resolved-IP sanity checks
    (no redirects followed, no private/loopback/link-local targets)."""
    from urllib.parse import urlparse
    parts = urlparse(url)
    if parts.scheme != "https" or parts.hostname not in ALLOWED_HOSTS:
        raise ValueError(f"URL outside the pinned allowlist: {url!r}")
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
    """Extract members whose names contain no traversal segments and
    whose resolved paths stay inside `dest`."""
    dest = dest.resolve()
    for member in tf.getmembers():
        name = member.name.replace("\\", "/")
        parts = Path(name).parts
        if (name.startswith(("/", "~")) or ".." in parts
                or Path(name).is_absolute()):
            raise ValueError(f"tar member rejected: {member.name!r}")
        target = (dest / member.name).resolve()
        if not str(target).startswith(str(dest)):
            raise ValueError(f"tar member escapes destination: {member.name!r}")
    tf.extractall(dest)  # members validated above


def ensure_checkout() -> Path:
    env = os.environ.get("BENCHPRESS_HOME")
    if env and Path(env).is_dir():
        return Path(env)
    local = REPO / "third_party" / "benchpress"
    if local.is_dir():
        return local
    local.parent.mkdir(exist_ok=True)
    print("downloading qiskit/benchpress main tarball ...", flush=True)
    data = _safe_fetch(BENCHPRESS_TARBALL)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tf:
        _safe_extract(tf, local.parent)
    extracted = local.parent / "benchpress-main"
    extracted.rename(local)
    return local


def install_gym(checkout: Path) -> None:
    dst = checkout / "benchpress" / "compactq_gym"
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(GYM_SRC, dst)
    sys.path.insert(0, str(PATCH.parent))
    import importlib.util
    spec = importlib.util.spec_from_file_location("apply_patches", PATCH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if mod.main(str(checkout)) != 0:
        raise SystemExit("patching benchpress failed")


def run_benchpress(checkout: Path, gym: str, sizes: list[str],
                   out_dir: Path, rounds: int) -> list[Path]:
    # repo first: the compactq under test is THIS checkout, not a
    # pip-installed copy; the checkout follows so `benchpress` imports
    sys.path.insert(0, str(REPO))
    sys.path.insert(1, str(checkout))
    written = []
    stamp = time.strftime("%Y%m%d-%H%M%S")
    gym_test = checkout / "benchpress" / f"{gym}_gym" / "abstract_transpile" / \
        "test_qasmbench.py"
    if not gym_test.is_dir() and not gym_test.is_file():
        raise SystemExit(f"gym test module missing: {gym_test}")
    for size in sizes:
        out_json = out_dir / f"benchpress_{gym}_{size}_{stamp}.json"
        workout = (checkout / "benchpress" / "workouts" /
                   "abstract_transpile" / "qasmbench.py")
        ini = checkout / "benchpress" / "pytest.ini"
        argv = [
            "-c", str(ini), str(workout), str(gym_test),
            f"--benchmark-json={out_json}",
            f"--benchmark-min-rounds={rounds}",
            "-k", f"TestWorkoutAbstractQasmBench{size.capitalize()} and all-to-all",
            "-p", "no:cacheprovider", "--no-header", "-q",
        ]
        print("RUN pytest:", " ".join(argv), flush=True)
        import pytest
        t0 = time.time()
        rc = pytest.main(argv)
        print(f"  -> exit {rc} in {time.time()-t0:.0f}s", flush=True)
        if out_json.exists():
            written.append(out_json)
        else:
            print(f"  WARNING: no record written for {gym}/{size}")
    return written


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gym", default="compactq", choices=("compactq", "qiskit"))
    ap.add_argument("--sizes", default="small,medium")
    ap.add_argument("--out-dir", default=str(REPO / "results" / "benchpress"))
    ap.add_argument("--rounds", type=int, default=1)
    args = ap.parse_args()
    sizes = [s.strip() for s in args.sizes.split(",") if s.strip()]
    if not sizes or any(s not in KNOWN_SIZES for s in sizes):
        print(f"invalid --sizes; expected from {sorted(KNOWN_SIZES)}")
        return 1
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    checkout = ensure_checkout()
    if args.gym == "compactq":
        install_gym(checkout)

    written = run_benchpress(checkout, args.gym, sizes, out_dir, args.rounds)

    meta_path = out_dir / "environment.json"
    meta = (json.loads(meta_path.read_text(encoding="utf-8"))
            if meta_path.exists() else {})
    meta.update({
        "generated": datetime.now(timezone.utc).isoformat(),
        "benchpress_checkout": str(checkout),
        "gym_files": "benchpress_integration/compactq_gym (MIT, compactq)",
        "protocol": ("Benchpress abstract-transpilation workout, QASMBench "
                     "suites, all-to-all topology, pytest-benchmark native "
                     "records; the timed region for compactq is exactly "
                     "optimize_search (its proof included)"),
    })
    meta["records"] = sorted(set(
        meta.get("records", []) + [p.name for p in written]))
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {meta_path} (+{len(written)} record files)")
    return 0 if written else 1


if __name__ == "__main__":
    sys.exit(main())
