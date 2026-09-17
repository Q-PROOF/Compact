"""SuppressionReport - the per-stage audit artifact.

Every closed suppression service returns numbers; this module returns the
*work*: which stages ran, what each one did to the circuit, whether the
result was proven (and at what proof level), and the measured suppression
factor when the built-in simulator executed the plan.

Proof-status vocabulary:
    exact-unitary   dense 4x4..256x256 fidelity proof (<= 8 qubits)
    tableau         Clifford-tableau equality (any size, Clifford circuits)
    randomized-K    randomized statevector agreement beyond the dense limit
    input-unchanged the stage's circuit is the untouched input (doubt path)
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class StageRecord:
    """One pipeline stage applied to one circuit."""
    name: str
    gates: int
    two_q: int
    depth: int
    proof: str
    note: str = ""


def zero_noise(num_qubits: int):
    """An explicitly silent NoiseModel (all calibrations zero) - used to
    derive the ideal outcome distribution for measured factors."""
    from .noise import NoiseModel
    m = NoiseModel(num_qubits)
    m.gate_infidelity = {"1q": 0.0, "2q": 0.0}
    for w in range(num_qubits):
        m.t1_us[w] = 1e12
        m.t2_us[w] = 1e12
        m.readout[w] = (0.0, 0.0)
    return m


class SuppressionReport:
    """Accumulates one StageRecord per pipeline stage plus execution
    metadata and (in simulator mode) the measured suppression factor."""

    def __init__(self, num_qubits: int):
        self.num_qubits = num_qubits
        self.stages: list = []
        self.meta: dict = {}
        self.measured: dict = {}

    # -- construction ------------------------------------------------------
    def add(self, name: str, circ, proof: str = "exact-unitary",
            note: str = "") -> None:
        self.stages.append(StageRecord(name, len(circ.ops),
                                       circ.two_qubit_count(), circ.depth(),
                                       proof, note))

    # -- measured results ---------------------------------------------------
    def set_measured(self, raw_prob: float, suppressed_prob: float,
                     ideal_support: list) -> None:
        factor = (suppressed_prob / raw_prob) if raw_prob > 0 else None
        self.measured = {"ideal_support": list(ideal_support),
                         "raw_ideal_probability": raw_prob,
                         "suppressed_ideal_probability": suppressed_prob,
                         "suppression_factor": (round(factor, 4)
                                                if factor is not None
                                                else None)}

    # -- serialization -------------------------------------------------------
    def to_dict(self) -> dict:
        d = {"num_qubits": self.num_qubits,
             "stages": [vars(s) for s in self.stages],
             "meta": dict(self.meta)}
        if self.measured:
            d["measured"] = dict(self.measured)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "SuppressionReport":
        r = cls(d["num_qubits"])
        r.stages = [StageRecord(**s) for s in d.get("stages", [])]
        r.meta = dict(d.get("meta", {}))
        r.measured = dict(d.get("measured", {}))
        return r

    # -- display --------------------------------------------------------------
    def __str__(self) -> str:
        rows = [f"{'stage':<14} {'gates':>6} {'2q':>4} {'depth':>6}  proof"]
        for s in self.stages:
            rows.append(f"{s.name:<14} {s.gates:>6} {s.two_q:>4} "
                        f"{s.depth:>6}  {s.proof}"
                        + (f" ({s.note})" if s.note else ""))
        if self.measured:
            m = self.measured
            f = m.get("suppression_factor")
            rows.append(f"\nmeasured: P(ideal) raw {m['raw_ideal_probability']:.4f}"
                        f" -> suppressed {m['suppressed_ideal_probability']:.4f}"
                        + (f"  (x{f})" if f else ""))
        for k, v in self.meta.items():
            rows.append(f"{k}: {v}")
        return "\n".join(rows)
