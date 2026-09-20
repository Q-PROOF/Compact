"""Interactive-latency benchmark: import time + small-circuit optimize
wall time (p50/p95), compared against Qiskit when installed.

Protocol: best-of-5 process spawns for import; 200 in-process optimize
calls on a 3-qubit / 12-gate circuit for wall-time percentiles.  The
proof is included in compactq's numbers (the CLI always proves).
"""
from __future__ import annotations

import json
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
# measure THIS repo, not a pip-installed copy
sys.path.insert(0, str(REPO))

IMPORT_PROBE = "import time; t0=time.perf_counter(); import compactq; " \
               "print(time.perf_counter()-t0)"
QISKIT_PROBE = "import time; t0=time.perf_counter(); import qiskit; " \
               "print(time.perf_counter()-t0)"


def best_import_time(probe: str, tries: int = 5) -> float:
    best = float("inf")
    for _ in range(tries):
        out = subprocess.run([sys.executable, "-c", probe],
                             capture_output=True, text=True,
                             cwd=str(REPO))
        if out.returncode == 0:
            best = min(best, float(out.stdout.strip()))
    return best


def main() -> int:
    import compactq
    from compactq import Circuit, Gate, optimize_search

    rng = __import__("random").Random(4)
    ops = []
    for _ in range(12):
        if rng.random() < 0.5:
            ops.append(Gate(rng.choice(["h", "t", "s"]), (), (rng.randrange(3),)))
        else:
            a, b = rng.sample(range(3), 2)
            ops.append(Gate("cx", (), (a, b)))
    circ = Circuit(3, ops)

    wall = []
    for _ in range(200):
        t0 = time.perf_counter()
        optimize_search(circ)
        wall.append((time.perf_counter() - t0) * 1000)

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": __import__("provenance").git_sha(),
        "protocol": ("import time: best-of-5 fresh-process spawns; "
                     "optimize wall: 200 calls on a fixed 3q/12-gate "
                     "circuit, proof included"),
        "compactq_import_ms": round(best_import_time(IMPORT_PROBE) * 1000, 2),
        "optimize_ms_p50": round(statistics.median(wall), 3),
        "optimize_ms_p95": round(sorted(wall)[int(0.95 * len(wall))], 3),
        "version": compactq.__version__,
    }
    try:
        doc["qiskit_import_ms"] = round(best_import_time(QISKIT_PROBE) * 1000, 2)
    except Exception:
        pass

    out = REPO / "results" / "latency.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(doc, indent=2))
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(REPO / "scripts"))
    raise SystemExit(main())
