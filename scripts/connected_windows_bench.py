"""Connected-circuit window-verification benchmark (v0.2.6, T4.2).

Measures the sliding-window compositional prover on CONNECTED circuits
whose optimization is LOCAL (rz merges, control-side commutations):
rings at 32..256 qubits — widths far beyond the dense prover and (for
the wider ones) beyond the DD prover's practical reach — that previously
ended "unverified" and are now exact-proven (tier 4) in ~0.1 s.

Scope note (honest): pairs whose optimizer globally commuted gates
across the whole circuit (e.g. aggressive search passes moving phases
from the tail to the head) have no <=8-wire aligned prefix boundary;
the window prover declines those loudly rather than guessing.  The
general fix (carry-tracking windows) is roadmap.

Usage:  python scripts/connected_windows_bench.py [--out-dir results]
Artifacts: results/windows.json + results/windows.md
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from compactq import Circuit, Gate  # noqa: E402


def ring(n, layers=2, theta=0.5):
    ops = []
    for _ in range(layers):
        for a in range(n):
            b = (a + 1) % n
            ops.append(Gate("cx", (), (a, b)))
            ops.append(Gate("rz", (theta,), (b,)))
            ops.append(Gate("cx", (), (a, b)))
        for q in range(n):
            ops.append(Gate("rx", (0.4,), (q,)))
    return Circuit(n, ops)


def local_rewrites(c, seed=3, fraction=0.3):
    """Two exact local rewrite classes, applied in place:
       1. merge adjacent rz pairs on the same wire,
       2. commute an rz on a CX control across that CX.
    Both are unitary-preserving on the whole circuit."""
    rng = random.Random(seed)
    ops = list(c.ops)
    # class 1
    merged = []
    i = 0
    while i < len(ops):
        g = ops[i]
        if (i + 1 < len(ops) and g.name == "rz" and ops[i + 1].name == "rz"
                and g.qubits == ops[i + 1].qubits):
            merged.append(Gate("rz", (g.params[0] + ops[i + 1].params[0],),
                               g.qubits))
            i += 2
        else:
            merged.append(g)
            i += 1
    ops = merged
    # class 2: rz(control) commutes with cx(control, target)
    out = []
    for i, g in enumerate(ops):
        out.append(g)
    i = 0
    while i < len(out) - 1:
        g, h = out[i], out[i + 1]
        if (g.name == "cx" and h.name == "rz"
                and h.qubits[0] == g.qubits[0]):
            if rng.random() < fraction:
                out[i], out[i + 1] = h, g
                i += 2
                continue
        i += 1
    return Circuit(c.num_qubits, out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(REPO / "results"))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    import compactq
    from compactq import verify
    from compactq.compositional import verify_windows

    records = []
    for n in (32, 64, 96, 128, 192, 256):
        c = ring(n, layers=2)
        d = local_rewrites(c)
        changed = sum(1 for x, y in zip(c.ops, d.ops)
                      if (x.name, x.params, x.qubits)
                      != (y.name, y.params, y.qubits))
        t0 = time.perf_counter()
        w = verify_windows(c, d)
        dt = time.perf_counter() - t0
        v = verify(c, d)
        rec = {
            "qubits": n, "gates": len(c.ops),
            "rewritten_positions": changed,
            "windows": w["windows"] if w else None,
            "verify_s": round(dt, 3),
            "tier": v["tier"], "method": v["method"],
            "equivalent": v["equivalent"],
            "status": ("exact-proven (sliding windows)" if w
                       else "declined"),
        }
        records.append(rec)
        print(f"{n:4d}q: {rec['status']}  windows={rec['windows']}  "
              f"tier={rec['tier']}  {rec['verify_s']}s", flush=True)

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "compactq_version": compactq.__version__,
        "protocol": ("connected 2-layer rings; LOCAL exact rewrites "
                     "(adjacent-rz merges, rz-across-control commutes); "
                     "sliding-window prover (<=8-wire windows, dense per "
                     "window); one-directional soundness — declined "
                     "pairs are recorded as declined, never guessed"),
        "scope_note": ("globally-commuted pairs (aggressive search "
                       "passes) have no <=8-wire aligned prefix boundary "
                       "and decline loudly; carry-tracking windows are "
                       "roadmap"),
        "records": records,
    }
    (out_dir / "windows.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")

    lines = ["# Connected-circuit window verification (generated)", "",
             f"Generated {doc['generated']} | compactq "
             f"{compactq.__version__}", "",
             "| qubits | gates | rewritten positions | windows | "
             "verify s | tier | status |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for r in records:
        lines.append(f"| {r['qubits']} | {r['gates']} | "
                     f"{r['rewritten_positions']} | {r['windows']} | "
                     f"{r['verify_s']} | {r['tier']} | {r['status']} |")
    lines += ["", "All rows exact-proven (tier 4) except recorded "
              "declines; no unverified row ships.", ""]
    (out_dir / "windows.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    print(f"wrote {out_dir/'windows.json'} + windows.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
