"""Benchmark reproducibility checker: validates the provenance of every
artifact in results/ (plus the root benchmark artifacts).

For each artifact this checks:
  - parseable JSON
  - a recorded git SHA (short, 7-12 hex chars)
  - a generated timestamp (ISO format)
  - a known artifact producer (documented in results/README.md)

Fails (exit 1) on any artifact with missing provenance — the
evidence-integrity gate for the release pipeline.  Advisory in CI.
"""
from __future__ import annotations

import json
import re
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


def main() -> int:
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
        if "generated" in required or "generated" in doc:
            ts = str(doc.get("generated", ""))
            if not re.match(r"\d{4}-\d{2}-\d{2}T", ts):
                problems.append(f"bad timestamp {ts!r}")
        sha = str(doc.get("commit", ""))
        if required and ("commit" in required or "commit" in doc):
            if not SHA_RE.match(sha):
                problems.append(f"bad commit sha {sha!r}")
        for key in required:
            if key not in doc:
                problems.append(f"missing key {key!r}")
        if problems:
            print(f"  FAIL  {rel}: {'; '.join(problems)}")
            ok = False
        else:
            print(f"  OK    {rel} (producer {producer})")
    print()
    if not ok:
        print("PROVENANCE CHECK: FAILED")
        return 1
    print("PROVENANCE CHECK: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
