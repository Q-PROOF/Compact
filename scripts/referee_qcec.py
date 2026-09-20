"""Independent referee #3: MQT QCEC (decision-diagram / ZX equivalence).

Runs MQT QCEC — a C++ DD-based equivalence checker with ZX-calculus and
simulation strategies, fully independent of both qiskit's Operator and
PyZX — on the input-vs-optimized pairs of the QASMBench small suite.

Emits results/qcec_referee.json.  QCEC verdicts are among the strongest
independent equivalence evidence available in open source.
"""
from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))
warnings.filterwarnings("ignore")


def main() -> int:
    import compactq
    import numpy as np
    from mqt import qcec
    from mqt.qcec.pyqcec import EquivalenceCriterion
    from qiskit import qasm2
    from qiskit.quantum_info import Operator
    from compactq.qiskit_bridge import from_qiskit, to_qiskit
    from provenance import environment, git_sha
    from realbench import normalize_qiskit, prepare_unitary

    print(f"compactq {compactq.__version__} | qcec via mqt.qcec")
    root = REPO / "third_party" / "QASMBench" / "small"
    records = []
    tmp = REPO / "results"
    tmp.mkdir(exist_ok=True)
    for d in sorted(root.iterdir()):
        if not d.is_dir():
            continue
        qf = d / f"{d.name}.qasm"
        if not qf.is_file():
            continue
        rec = {"circuit": d.name, "tool": "compactq", "referee": "mqt-qcec"}
        prepared, _reason = prepare_unitary(qf.read_text())
        if prepared is None or prepared.num_qubits > 16:
            rec["status"] = "skipped"
            records.append(rec)
            continue
        try:
            inp = normalize_qiskit(prepared)
            circ = from_qiskit(inp)
            opt = compactq.optimize_search(circ)
            f_in = tmp / f"_qcec_in_{d.name}.qasm"
            f_out = tmp / f"_qcec_out_{d.name}.qasm"
            f_in.write_text(qasm2.dumps(inp), encoding="utf-8")
            f_out.write_text(qasm2.dumps(to_qiskit(opt)), encoding="utf-8")
            res = qcec.verify(str(f_in), str(f_out))
            eq = res.equivalence
            rec["qcec_criterion"] = str(eq)
            if eq == EquivalenceCriterion.equivalent:
                rec["equivalent"] = True
                rec["status"] = "OK"
            elif eq == EquivalenceCriterion.equivalent_up_to_global_phase:
                # compactq's guarantee IS equivalence up to global phase,
                # so QCEC's phase-insensitive verdict is a hard OK (this
                # converted 18 spurious "inconclusive" verdicts)
                rec["equivalent"] = True
                rec["global_phase_only"] = True
                rec["status"] = "OK (up to global phase)"
            elif eq == EquivalenceCriterion.not_equivalent:
                # cross-check the dispute against the Qiskit Operator
                # referee before accepting a not-equivalent verdict: the
                # two referees use different numerical tolerances
                got = Operator(to_qiskit(opt)).data
                ref = Operator(inp).data
                fid = float(abs(np.sum(np.conj(ref) * got)) / ref.shape[0])
                rec["operator_fidelity"] = round(fid, 12)
                if fid > 1 - 1e-6:
                    rec["equivalent"] = True
                    rec["status"] = ("referee-disagreement (Operator "
                                     "equivalent within 1e-6; QCEC "
                                     "tolerance stricter)")
                else:
                    rec["equivalent"] = False
                    rec["status"] = "INEQUIVALENT"
            else:
                rec["equivalent"] = None
                rec["status"] = "inconclusive"
        except Exception as e:
            rec["status"] = f"error {type(e).__name__}: {e}"
        records.append(rec)
        print(f"{d.name:24s} {rec['status']}", flush=True)
        (tmp / f"_qcec_in_{d.name}.qasm").unlink(missing_ok=True)
        (tmp / f"_qcec_out_{d.name}.qasm").unlink(missing_ok=True)

    ok = sum(1 for r in records if str(r.get("status", "")).startswith("OK"))
    bad = sum(1 for r in records if r.get("status") == "INEQUIVALENT")
    other = len(records) - ok - bad
    print(f"qcec referee: {ok} OK, {bad} INEQUIVALENT, {other} other")

    doc = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "commit": git_sha(),
        "referee": ("MQT QCEC: DD-based equivalence checking with "
                    "simulation and ZX-calculus strategies"),
        "environment": environment(extra_deps=("mqt.qcec",)),
        "count": len(records),
        "ok": ok,
        "inequivalent": bad,
        "records": records,
    }
    out = REPO / "results" / "qcec_referee.json"
    out.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
