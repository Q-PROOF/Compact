"""Same-input, same-gate-set reproduction harness.

Closes the benchmark-fairness gap identified in the 2026-09-18 external
audit: every tool receives THE SAME input circuit (generated here in
qiskit as the common interchange), every output is LOWERED TO THE SAME
GATE SET (u3 + cx) before counting, counts are taken at two levels
(CX-native 2q, and CP-native where cp counts as one 2-qubit gate), and
every output is refereed with qiskit's Operator.

The suite is PINNED by results/MANIFEST.json: per-circuit QASM digests
plus referee/tool versions.  Any drift (a qiskit version change that
alters a library circuit, a changed generator) fails the run loudly
instead of silently producing numbers a stranger cannot reproduce.
Regenerate the manifest deliberately with --update-manifest.

Usage:

    python scripts/repro_harness.py [--tools compact,qiskit,pytket,cirq]
                                    [--out-dir results]
                                    [--update-manifest]

Artifacts: results/repro.json + results/repro.csv + results/repro.md
Pinning:   results/MANIFEST.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
warnings.filterwarnings("ignore")

# tool outputs are lowered before counting: this harness is the referee,
# the in-product proof is not what is being measured here
_OFF = False

MANIFEST_VERSION = 1
REF_TOL = 1e-6

import numpy as np  # noqa: E402
from qiskit import QuantumCircuit, qasm2, transpile  # noqa: E402
from qiskit.quantum_info import Operator  # noqa: E402

from provenance import environment, git_sha  # noqa: E402

BASIS = ["u3", "cx"]


def _qasm_digest(qc: QuantumCircuit) -> str:
    return hashlib.sha256(qasm2.dumps(qc).encode("utf-8")).hexdigest()


def build_suite() -> dict[str, QuantumCircuit]:
    """Deterministic suite mirroring the external audit's families."""
    from qiskit.circuit.library import EfficientSU2, QFT
    suite: dict[str, QuantumCircuit] = {}

    for n in (4, 6, 8):
        suite[f"qft_{n}"] = QFT(n, do_swaps=True).decompose()

    for n in (4, 6, 8):  # QAOA-like: cost + mixer layers
        qc = QuantumCircuit(n)
        for _ in range(2):
            for q in range(n - 1):
                qc.rzz(0.8, q, q + 1)
            for q in range(n):
                qc.rx(1.1, q)
        suite[f"qaoa_{n}"] = qc

    for n in (4, 6):  # QPE-like: counting register + phase ladder + invQFT
        import math
        qc = QuantumCircuit(n)
        counting = max(2, n // 2)
        qc.h(range(counting))
        qc.x(n - 1)
        for j in range(counting):
            for _ in range(2 ** j):
                qc.cp(math.pi / 4, j, n - 1)
        for j in range(counting // 2):
            qc.swap(j, counting - 1 - j)
        for j in range(counting):
            for k in range(j):
                qc.cp(-math.pi / 2 ** (j - k), k, j)
            qc.h(j)
        suite[f"qpe_{n}"] = qc

    for n in (3, 5):  # Grover: oracle phase flip + diffuser on n qubits
        qc = QuantumCircuit(n)
        for q in range(n):
            qc.h(q)
        qc.compose(grover_iter(n), inplace=True)
        qc.compose(grover_iter(n), inplace=True)
        suite[f"grover_{n}"] = qc

    for n in (4, 6):  # Heisenberg/Trotter: XX+YY+ZZ evolution, 2 layers
        qc = QuantumCircuit(n)
        for _ in range(2):
            for q in range(n - 1):
                qc.rxx(0.6, q, q + 1)
                qc.ryy(0.4, q, q + 1)
                qc.rzz(0.5, q, q + 1)
        suite[f"heisenberg_{n}"] = qc

    for n in (5,):  # Clifford and Clifford+T stress
        import random
        rng = random.Random(n * 17)
        qc = QuantumCircuit(n)
        for _ in range(4 * n):
            q = rng.randrange(n)
            g = rng.choice(("h", "s", "x", "z"))
            getattr(qc, g)(q)
            if rng.random() < 0.5 and n >= 2:
                a = rng.randrange(n)
                b = (a + rng.randrange(1, n)) % n
                qc.cx(a, b)
        suite[f"clifford_{n}"] = qc
        qc2 = QuantumCircuit(n)
        for _ in range(2 * n):
            q = rng.randrange(n)
            qc.t(q)
            qc.h(q)
            qc.s(q)
            if n >= 2:
                a = rng.randrange(n)
                b = (a + rng.randrange(1, n)) % n
                qc.t(b)
                qc.cx(a, b)
        suite[f"cliffordt_{n}"] = qc2

    for n in (4, 6):  # ripple-carry style adder via Cuccaro-ish MCX ladder
        qc = QuantumCircuit(n)
        for q in range(n - 1):
            qc.ccx(q, q + 1, (q + 2) % n) if q + 2 < n else qc.cx(q, q + 1)
        qc.cx(n - 1, 0)
        suite[f"adder_{n}"] = qc

    for n in (4, 6):  # VQE hardware-efficient ansatz (parameters bound)
        import numpy as np
        ansatz = EfficientSU2(n, reps=1, entanglement="linear")
        qc = ansatz.decompose()
        qc.assign_parameters(np.linspace(0.3, 1.1, qc.num_parameters),
                             inplace=True)
        suite[f"vqe_{n}"] = qc

    for n in (4, 6):  # redundancy stress: repeated self-inverse cascades
        qc = QuantumCircuit(n)
        for rep in range(4):
            for q in range(n):
                qc.t(q)
                qc.tdg(q)
            for q in range(n - 1):
                qc.cx(q, q + 1)
            for q in reversed(range(n - 1)):
                qc.cx(q, q + 1)
        suite[f"redundant_{n}"] = qc

    return suite


def grover_iter(n: int) -> QuantumCircuit:
    qc = QuantumCircuit(n)
    qc.h(n - 1)
    for q in range(n):  # oracle: phase flip on |11..1>
        qc.x(q)
    qc.h(n - 2)
    qc.mcx(list(range(n - 2)) + [n - 1], n - 2) if n > 2 else qc.cz(0, 1)
    qc.h(n - 2)
    for q in range(n):
        qc.x(q)
    qc.h(n - 1)
    # diffuser
    for q in range(n):
        qc.h(q)
        qc.x(q)
    qc.h(n - 1)
    if n > 2:
        qc.mcx(list(range(n - 1)), n - 1)
    elif n == 2:
        qc.cz(0, 1)
    qc.h(n - 1)
    for q in range(n):
        qc.x(q)
        qc.h(q)
    return qc


def lower_2q(qc: QuantumCircuit) -> QuantumCircuit:
    """Lower to the common gate set u3+cx, no further optimization."""
    return transpile(qc, basis_gates=BASIS, optimization_level=0)


def count_2q(qc: QuantumCircuit) -> int:
    return sum(1 for i in qc.data if len(i.qubits) == 2)


def count_cp_native(qc: QuantumCircuit) -> int:
    """CP-native view: cx + cp both count as one 2-qubit gate."""
    return sum(1 for i in qc.data if len(i.qubits) == 2
               and i.operation.name in ("cx", "cp"))


def compact_run(inp: QuantumCircuit):
    import compactq
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    circ = from_qiskit(inp)
    t0 = time.perf_counter()
    out = compactq.optimize_search(circ, verify=_OFF)
    dt = time.perf_counter() - t0
    return to_qiskit(out), dt


def qiskit_run(inp: QuantumCircuit):
    t0 = time.perf_counter()
    out = transpile(inp, optimization_level=3)
    dt = time.perf_counter() - t0
    return out, dt


def pytket_run(inp: QuantumCircuit):
    from pytket.extensions.qiskit import qiskit_to_tk, tk_to_qiskit
    from pytket.passes import FullPeepholeOptimise
    t0 = time.perf_counter()
    tk = qiskit_to_tk(inp)
    FullPeepholeOptimise().apply(tk)
    out = tk_to_qiskit(tk)
    dt = time.perf_counter() - t0
    return out, dt


def cirq_run(inp: QuantumCircuit):
    import cirq
    from cirq.contrib.qasm_import import circuit_from_qasm
    t0 = time.perf_counter()
    cq = circuit_from_qasm(qasm2.dumps(inp))
    optimized = cirq.optimize_for_target_gateset(
        cq, gateset=cirq.CZTargetGateset())
    dt = time.perf_counter() - t0
    out = QuantumCircuit.from_qasm_str(cirq.qasm(optimized))
    return out, dt


TOOL_RUNNERS = {"compact": compact_run, "qiskit": qiskit_run,
                "pytket": pytket_run, "cirq": cirq_run}


def pyzx_run(inp: QuantumCircuit):
    import pyzx
    t0 = time.perf_counter()
    circ = pyzx.Circuit.from_qasm(qasm2.dumps(inp))
    graph = circ.to_graph()
    pyzx.full_reduce(graph, quiet=True)
    out_qasm = pyzx.extract_circuit(graph.copy()).to_qasm()
    dt = time.perf_counter() - t0
    return QuantumCircuit.from_qasm_str(out_qasm), dt


def bqskit_run(inp: QuantumCircuit):
    import bqskit
    from bqskit import compile as bqskit_compile
    t0 = time.perf_counter()
    bqs = bqskit.Circuit.from_qasm(qasm2.dumps(inp))
    compiled = bqskit_compile(bqs)
    out = compiled.to_qiskit()
    dt = time.perf_counter() - t0
    return out, dt


TOOL_RUNNERS["pyzx"] = pyzx_run
TOOL_RUNNERS["bqskit"] = bqskit_run


def _manifest(out_dir: Path, suite: dict, tools: list) -> dict:
    """Build the manifest dict describing the suite + referee pins."""
    import qiskit
    circuits = {}
    for name, inp in suite.items():
        circuits[name] = {"qubits": inp.num_qubits,
                          "depth": inp.depth(),
                          "qasm_sha256": _qasm_digest(inp),
                          "lowered_sha256": _qasm_digest(lower_2q(inp))}
    return {
        "manifest_version": MANIFEST_VERSION,
        "protocol": ("same input for every tool; every output lowered to "
                     "u3+cx before counting; every output refereed by "
                     "qiskit.quantum_info.Operator, fidelity tolerance "
                     f"{REF_TOL}"),
        "referee": {"tool": "qiskit.quantum_info.Operator",
                    "qiskit_version_pinned": qiskit.__version__,
                    "tolerance": REF_TOL},
        "tools": tools,
        "suite_generator": "scripts/repro_harness.py::build_suite",
        "circuits": circuits,
        "outputs": ["repro.json", "repro.csv", "repro.md"],
    }


def _check_manifest(out_dir: Path, manifest: dict, update: bool) -> int:
    """Compare the freshly built suite against the committed pin.
    Returns 0 when compatible, 1 on mismatch (unless updating)."""
    path = out_dir / "MANIFEST.json"
    if not path.exists():
        path.write_text(json.dumps(manifest, indent=2,
                                   sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {path} ({len(manifest['circuits'])} circuits pinned)")
        return 0
    old = json.loads(path.read_text(encoding="utf-8"))
    problems = []
    if old.get("manifest_version") != MANIFEST_VERSION:
        problems.append(f"manifest_version {old.get('manifest_version')} "
                        f"!= {MANIFEST_VERSION}")
    old_ref = (old.get("referee") or {}).get("qiskit_version_pinned")
    new_ref = manifest["referee"]["qiskit_version_pinned"]
    if old_ref != new_ref:
        problems.append(f"referee qiskit {old_ref} -> {new_ref}")
    for name, new in manifest["circuits"].items():
        prev = old.get("circuits", {}).get(name)
        if prev is None:
            problems.append(f"{name}: newly added to the suite")
        elif prev.get("qasm_sha256") != new["qasm_sha256"] or \
                prev.get("lowered_sha256") != new["lowered_sha256"]:
            problems.append(f"{name}: circuit digest changed "
                            f"(generator/library drift)")
    for name in sorted(set(old.get("circuits", {})) -
                       set(manifest["circuits"])):
        problems.append(f"{name}: removed from the suite")
    if problems:
        if update:
            path.write_text(json.dumps(manifest, indent=2,
                                       sort_keys=True) + "\n",
                            encoding="utf-8")
            print(f"MANIFEST UPDATED ({len(problems)} changes):")
        else:
            print("MANIFEST MISMATCH — numbers would not be comparable "
                  "with the committed pin:")
        for p in problems:
            print(f"  - {p}")
        if update:
            print(f"rewrote {path}")
            return 0
        print("fix the environment or re-pin deliberately with "
              "--update-manifest")
        return 1
    print(f"manifest OK: {len(manifest['circuits'])} circuits, referee "
          f"qiskit {new_ref}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="same-input reproduction harness")
    ap.add_argument("--tools", default="compact,qiskit,pytket,cirq")
    ap.add_argument("--out-dir", default=str(REPO / "results"))
    ap.add_argument("--update-manifest", action="store_true",
                    help="re-pin results/MANIFEST.json to this environment")
    args = ap.parse_args()
    tools = [t.strip() for t in args.tools.split(",") if t.strip()]
    for t in tools:
        if t not in TOOL_RUNNERS:
            print(f"unknown tool {t!r}; expected {sorted(TOOL_RUNNERS)}")
            return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    suite = build_suite()
    rc = _check_manifest(out_dir, _manifest(out_dir, suite, tools),
                         update=args.update_manifest)
    if rc:
        return rc

    print(f"suite: {len(suite)} circuits; tools: {tools}")
    records = []
    for name, inp in suite.items():
        inp_l = lower_2q(inp)
        ref = Operator(inp_l).data
        in_2q = count_2q(inp_l)
        row = {"circuit": name, "qubits": inp.num_qubits,
               "input_2q_cx_lowered": in_2q,
               "input_2q_cp_native": count_cp_native(inp)}
        for tool in tools:
            try:
                out, dt = TOOL_RUNNERS[tool](inp)
                low = lower_2q(out)
                got = Operator(low).data
                fid = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                row[tool] = {"gates": sum(1 for _ in low.data),
                             "two_qubit": count_2q(low),
                             "depth": low.depth(),
                             "output_fidelity": round(fid, 12),
                             "status": "OK" if fid > 1 - 1e-6
                             else "INEQUIVALENT",
                             "wall_ms": round(dt * 1000)}
            except Exception as e:
                row[tool] = {"status": f"error {type(e).__name__}: {e}"}
        records.append(row)
        line = f"{name:16s}"
        for tool in tools:
            r = row.get(tool, {})
            line += f" | {tool}: {r.get('two_qubit', '?'):>3} 2q"
        print(line, flush=True)

    summary = {}
    for tool in tools:
        two = [r[tool]["two_qubit"] for r in records
               if isinstance(r.get(tool), dict)
               and r[tool].get("status") == "OK"]
        cuts = [(r["input_2q_cx_lowered"] - r[tool]["two_qubit"])
                / r["input_2q_cx_lowered"] * 100
                for r in records if isinstance(r.get(tool), dict)
                and r[tool].get("status") == "OK"
                and r["input_2q_cx_lowered"] > 0]
        med = statistics.median([r[tool]["wall_ms"] for r in records
                                 if isinstance(r.get(tool), dict)
                                 and r[tool].get("wall_ms") is not None])
        summary[tool] = {
            "mean_2q_cut_pct": round(statistics.mean(cuts), 1) if cuts else None,
            "median_wall_ms": round(med),
        }
    w = t = l = 0
    for r in records:
        if not (isinstance(r.get("compact"), dict)
                and isinstance(r.get("qiskit"), dict)):
            continue
        if r["compact"].get("status") != "OK" or r["qiskit"].get("status") != "OK":
            continue
        if r["compact"]["two_qubit"] < r["qiskit"]["two_qubit"]:
            w += 1
        elif r["compact"]["two_qubit"] == r["qiskit"]["two_qubit"]:
            t += 1
        else:
            l += 1
    summary["compact_vs_qiskit_2q_WTL"] = [w, t, l]

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "method": ("same input for every tool; every output lowered to "
                   "u3+cx before counting; every output refereed by "
                   "qiskit.quantum_info.Operator"),
        "suite": sorted(suite.keys()),
        "summary": summary,
        "records": records,
        "environment": environment(extra_deps=("pytket", "cirq", "bqskit")),
    }
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)
    (out_dir / "repro.json").write_text(
        json.dumps(doc, indent=2) + "\n", encoding="utf-8")

    with (out_dir / "repro.csv").open("w", newline="", encoding="utf-8") as fh:
        cols = ["circuit", "qubits", "input_2q_cx_lowered",
                "input_2q_cp_native"]
        for t in tools:
            cols += [f"{t}_gates", f"{t}_2q", f"{t}_depth", f"{t}_fidelity",
                     f"{t}_status", f"{t}_wall_ms"]
        wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in records:
            flat = {k: v for k, v in r.items() if not isinstance(v, dict)}
            for t in tools:
                tr = r.get(t) or {}
                if isinstance(tr, dict):
                    flat[f"{t}_gates"] = tr.get("gates")
                    flat[f"{t}_2q"] = tr.get("two_qubit")
                    flat[f"{t}_depth"] = tr.get("depth")
                    flat[f"{t}_fidelity"] = tr.get("output_fidelity")
                    flat[f"{t}_status"] = tr.get("status")
                    flat[f"{t}_wall_ms"] = tr.get("wall_ms")
            wr.writerow(flat)

    lines = ["# Independent reproduction harness (generated)", "",
             f"Generated {doc['generated']} | commit `{doc['commit']}` | "
             f"tools: {', '.join(tools)}", "",
             "| circuit | qubits | in 2q (u3+cx) | " +
             " | ".join(f"{t} 2q" for t in tools) + " | " +
             " | ".join(f"{t} status" for t in tools) + " |",
             "|---|---:|---:| " + "|".join("---:" for _ in tools) + "| " +
             "|".join("---" for _ in tools) + "|"]
    for r in records:
        cells = " | ".join(
            str(r.get(t, {}).get("two_qubit", "?")) if isinstance(r.get(t), dict)
            else "?" for t in tools)
        stats = " | ".join(
            str(r.get(t, {}).get("status", "?")) if isinstance(r.get(t), dict)
            else "?" for t in tools)
        lines.append(f"| {r['circuit']} | {r['qubits']} | "
                     f"{r['input_2q_cx_lowered']} | {cells} | {stats} |")
    lines += ["", "## Summary", ""]
    for tool, s in summary.items():
        if isinstance(s, dict):
            lines.append(f"- **{tool}**: mean 2q cut "
                         f"{s['mean_2q_cut_pct']}%, median wall "
                         f"{s['median_wall_ms']} ms")
        else:
            lines.append(f"- compact vs qiskit 2q W/T/L: {s}")
    (out_dir / "repro.md").write_text("\n".join(lines) + "\n",
                                      encoding="utf-8")
    print(f"wrote {out_dir / 'repro.json'} + repro.csv + repro.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
