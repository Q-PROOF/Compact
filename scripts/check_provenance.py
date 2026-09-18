"""Benchmark reproducibility checker — two modes.

Default (schema mode): validates the provenance of every artifact in
results/ (plus the root benchmark artifacts):
  - parseable JSON
  - a recorded git SHA (short, 7-12 hex chars)
  - a generated timestamp (ISO format)
  - a known artifact producer (documented in results/README.md)

`--strict-current-sha`: additionally requires every present artifact's
recorded commit to equal the CURRENT repository HEAD.  Intended for use
immediately after a regeneration pass (e.g. in the release pipeline) —
artifacts legitimately lag HEAD between releases, so CI runs the default
mode as a BLOCKING gate while strict mode is opt-in.

Fails (exit 1) on any artifact with missing provenance (default) or any
SHA mismatch (strict).
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SHA_RE = re.compile(r"^[0-9a-f]{7,12}$")

# artifact relative path -> (producer script, required keys)
ARTIFACTS = {
    "bench_results.json": ("scripts/bench_json.py",
                           ("generated", "commit", "package_version")),
    "head_to_head_results.json": ("scripts/head_to_head.py", ()),
    "suppress_results.json": ("scripts/suppress_bench.py", ()),
    "dd_sequence_results.json": ("scripts/suppress_bench.py", ()),
    "scale_results.json": ("scripts/scale_bench.py", ()),
    "results/bqskit.json": ("scripts/bqskit_bench.py",
                            ("generated", "commit", "records")),
    "results/pytket_referee.json": ("scripts/referee_pytket.py",
                                    ("generated", "commit", "records")),
    "results/scalability.json": ("scripts/scale_gauntlet.py",
                                 ("generated", "commit", "compactq_version",
                                  "environment", "ceilings", "records")),
}


def current_head_sha() -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                          text=True, cwd=str(REPO)).stdout.strip()[:12]


def main() -> int:
    ap = argparse.ArgumentParser(description="artifact provenance gate")
    ap.add_argument("--strict-current-sha", action="store_true",
                    help="also require every present artifact's recorded "
                         "commit to equal the current HEAD (use right after "
                         "a regeneration pass)")
    args = ap.parse_args()
    head = current_head_sha() if args.strict_current_sha else None

    ok = True
    for rel, (producer, required) in ARTIFACTS.items():
        path = REPO / rel
        if not path.is_file():
            print(f"  SKIP  {rel} (not present)")
            continue
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"  FAIL  {rel}: unparseable JSON ({e})")
            ok = False
            continue
        problems = []
        ts = str(doc.get("generated", ""))
        if not re.match(r"\d{4}-\d{2}-\d{2}T", ts):
            problems.append(f"bad timestamp {ts!r}")
        sha = str(doc.get("commit", ""))
        if not SHA_RE.match(sha):
            problems.append(f"bad commit sha {sha!r}")
        elif args.strict_current_sha and sha != head:
            problems.append(f"artifact sha {sha!r} != current HEAD {head!r}")
        for key in required:
            if key not in doc:
                problems.append(f"missing key {key!r}")
        if problems:
            print(f"  FAIL  {rel}: {'; '.join(problems)}")
            ok = False
        else:
            extra = f" (sha matches HEAD)" if args.strict_current_sha else ""
            print(f"  OK    {rel} (producer {producer}){extra}")
    print()
    if not ok:
        print("PROVENANCE CHECK: FAILED")
        return 1
    print("PROVENANCE CHECK: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
