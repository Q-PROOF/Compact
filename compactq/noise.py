"""Hardware noise models for the error-suppression pipeline.

A NoiseModel carries everything the suppression passes need, with sane
superconducting defaults so every pass also runs standalone:

    t1_us / t2_us       per-wire relaxation times (microseconds)
    gate_infidelity     {(gate_name, wires): infidelity}  (1 - F)
    readout             per-wire (P(1|0), P(0|1)) confusion probabilities
    drift_rate          per-wire rms quasi-static Z rate (rad/ns) - the
                        low-frequency noise dynamical decoupling refocuses
    durations_ns        {gate_name: duration in ns} for scheduling/DD

`from_qiskit_backend` ingests live calibration data when qiskit-ibm-runtime
is available; `from_dict` covers every other vendor (the fields are plain
data).  The model is intentionally read-only evidence: passes consume it,
never mutate it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

_DEFAULT_DURATIONS_NS = {"1q": 50.0, "cx": 300.0, "cz": 300.0,
                         "swap": 900.0, "ecr": 300.0, "iswap": 300.0,
                         "cp": 300.0}


@dataclass
class NoiseModel:
    num_qubits: int
    t1_us: dict = field(default_factory=dict)      # wire -> T1 (us)
    t2_us: dict = field(default_factory=dict)      # wire -> T2 (us)
    gate_infidelity: dict = field(default_factory=dict)   # (name, wires) -> 1-F
    readout: dict = field(default_factory=dict)    # wire -> (p10, p01)
    drift_rate: dict = field(default_factory=dict)  # wire -> rms quasi-static Z rate (rad/ns)
    durations_ns: dict = field(default_factory=lambda: dict(_DEFAULT_DURATIONS_NS))
    coherent_fraction: float = 1.0   # share of calibrated gate error treated as coherent

    # -- convenience -------------------------------------------------------
    def t1(self, w: int) -> float:
        return self.t1_us.get(w, 200.0)

    def t2(self, w: int) -> float:
        return min(self.t2_us.get(w, 100.0), 2 * self.t1(w))

    def drift(self, w: int) -> float:
        """RMS quasi-static Z-rotation rate (rad/ns) on wire w (0 if unknown).
        This is the low-frequency noise dynamical decoupling refocuses."""
        return self.drift_rate.get(w, 0.0)

    def duration_of(self, g) -> float:
        """Duration (ns) of a Gate: explicit entry, else width-based default."""
        key = (g.name, g.qubits)
        if key in self.durations_ns:
            return float(self.durations_ns[key])
        if g.name in self.durations_ns:
            return float(self.durations_ns[g.name])
        d = self.durations_ns
        if len(g.qubits) == 1:
            return float(d.get("1q", 50.0))
        return float(d.get(g.name, d.get("2q", 300.0)))

    def gate_error(self, g) -> float:
        """Infidelity of a Gate: per-(name, wires) entry, else per-name,
        else the width default ("1q" / "2q") - mirroring duration_of.  For
        two-qubit gates the REVERSED pair is honored too (edges are
        undirected for pricing; direction matters only to the compiler)."""
        key = (g.name, tuple(g.qubits))
        if key in self.gate_infidelity:
            return self.gate_infidelity[key]
        if g.name in self.gate_infidelity:
            return self.gate_infidelity[g.name]
        if len(g.qubits) == 2:
            rev = self.gate_infidelity.get((g.name, tuple(reversed(g.qubits))))
            if rev is not None:
                return rev
            return float(self.gate_infidelity.get("2q", 0.008))
        return float(self.gate_infidelity.get("1q", 0.0005))

    def readout_error(self, w: int):
        """(P(1|0), P(0|1)) for wire w; defaults to (1%, 1%)."""
        p = self.readout.get(w, (0.01, 0.01))
        return (float(p[0]), float(p[1]))

    # -- ingestion ----------------------------------------------------------
    @classmethod
    def from_dict(cls, num_qubits: int, data: dict) -> "NoiseModel":
        m = cls(num_qubits)
        m.t1_us = {int(k): float(v) for k, v in (data.get("t1_us") or {}).items()}
        m.t2_us = {int(k): float(v) for k, v in (data.get("t2_us") or {}).items()}
        m.readout = {int(k): (float(v[0]), float(v[1]))
                     for k, v in (data.get("readout") or {}).items()}
        m.drift_rate = {int(k): float(v)
                        for k, v in (data.get("drift_rate") or {}).items()}
        m.coherent_fraction = float(data.get("coherent_fraction", 1.0))
        for k, v in (data.get("gate_infidelity") or {}).items():
            if isinstance(k, tuple):
                m.gate_infidelity[(k[0], tuple(k[1]))] = float(v)
            else:
                m.gate_infidelity[k] = float(v)
        m.durations_ns.update({k: float(v)
                               for k, v in (data.get("durations_ns") or {}).items()})
        return m

    @classmethod
    def from_qiskit_backend(cls, backend) -> "NoiseModel":
        """Ingest calibration data from a qiskit-ibm-runtime BackendV2
        (or a fake provider backend).  Falls back to defaults per field."""
        n = backend.num_qubits
        m = cls(n)
        props = getattr(backend, "properties", None)
        if callable(props):
            try:
                props = props()
            except Exception:
                props = None
        if props is not None:
            try:
                for i, q in enumerate(props.qubits):
                    vals = {item.name: item.value for item in q
                            if hasattr(item, "name")}
                    if "T1" in vals:
                        m.t1_us[i] = float(vals["T1"]) / 1000.0  # ns -> us
                    if "T2" in vals:
                        m.t2_us[i] = float(vals["T2"]) / 1000.0
                    if "readout_error" in vals:
                        e = float(vals["readout_error"])
                        m.readout[i] = (e, e)
                for g in props.gates:
                    infid = None
                    for item in g.parameters:
                        if getattr(item, "name", "") == "gate_error":
                            infid = float(item.value)
                    wires = tuple(g.qubits)
                    for name in (g.gate, getattr(g, "name", g.gate)):
                        if infid is not None:
                            m.gate_infidelity[(name, wires)] = infid
            except Exception:
                pass  # keep defaults for whatever failed to parse
        target = getattr(backend, "target", None)
        duration_map = getattr(target, "duration", None) \
            if target is not None else None
        if duration_map:
            try:
                for name, prop in duration_map.items():
                    if prop and prop[0] is not None:
                        m.durations_ns[name] = float(prop[0])  # already ns (dt*0.5)
            except Exception:
                pass  # keep defaults for whatever failed to parse
        return m


def default_model(num_qubits: int) -> NoiseModel:
    """A representative superconducting device: T1 150us, T2 80us, 1% readout,
    0.5% two-qubit / 0.05% one-qubit infidelity."""
    m = NoiseModel(num_qubits)
    for w in range(num_qubits):
        m.t1_us[w] = 150.0 + 20.0 * ((w * 7) % 5 - 2)
        m.t2_us[w] = 80.0 + 10.0 * ((w * 3) % 4 - 1)
        m.readout[w] = (0.01 + 0.004 * (w % 3), 0.012 + 0.003 * (w % 2))
    m.gate_infidelity = {"1q": 0.0005, "cx": 0.008, "cz": 0.008,
                         "swap": 0.024, "ecr": 0.008, "iswap": 0.008}
    return m
