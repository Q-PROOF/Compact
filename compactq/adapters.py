"""Thin execution adapters - the `execute(...)` box of the pipeline.

Turn a hardware handle into the pair the suppression pipeline consumes:
a `run_fn(circuit, *, seed, shots) -> counts` callable plus a NoiseModel
ingested from the same device's calibrations.

    from compactq import suppress_execute
    from compactq.adapters import qiskit_runtime

    rt = qiskit_runtime(backend)                 # BackendV2, incl. fakes
    result = suppress_execute(circ, rt.noise_model, rt.run_fn)

Every import beyond the standard library is lazy: the CORE stays
zero-dependency, and adapters raise a clear RuntimeError when the
optional stack (qiskit-ibm-runtime, amazon-braket-sdk) is missing.

Counts convention: hardware APIs return most-significant-bit-first
bitstrings; the adapters reverse them into compactq's little-endian
convention (char i of the string is wire i).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


@dataclass
class RunTarget:
    """A ready-to-use execution target: callable + calibrations."""
    run_fn: Callable
    noise_model: object
    name: str = "device"
    meta: dict = field(default_factory=dict)


def _reverse_counts(raw: dict, n: int) -> dict:
    """qiskit/braket count strings are MSB-first; compactq wants char i =
    wire i.  Truncate to n wires and reverse."""
    out = {}
    for s, c in raw.items():
        clean = s.replace(" ", "")
        if len(clean) < n:
            continue
        out[clean[::-1][:n]] = out.get(clean[::-1][:n], 0) + c
    return out


def qiskit_runtime(backend) -> RunTarget:
    """Wrap a qiskit BackendV2 (real via qiskit-ibm-runtime, or a
    fake-provider backend) into a RunTarget.  Uses SamplerV2; measurements
    are appended to a copy of each submitted circuit."""
    from .noise import NoiseModel
    noise = NoiseModel.from_qiskit_backend(backend)
    name = getattr(backend, "name", None) or str(getattr(backend,
                                                         "backend_name", ""))
    try:
        from qiskit_ibm_runtime import SamplerV2
        sampler = SamplerV2(mode=backend)
    except Exception as e:
        raise RuntimeError(
            f"qiskit-ibm-runtime SamplerV2 unavailable for {name}: {e}. "
            "Install qiskit-ibm-runtime or pass run_fn manually.") from e

    def run_fn(circ, *, seed: int = 0, shots: int = 1024) -> dict:
        from .qiskit_bridge import to_qiskit
        qc = to_qiskit(circ)
        qc.measure_all()
        try:
            job = sampler.run([qc], shots=shots)
        except TypeError:
            job = sampler.run([qc])          # some backends fix shots elsewhere
        result = job.result()
        data = result[0].data
        raw = data.meas.get_counts() if hasattr(data, "meas") \
            else data.get_counts()
        return _reverse_counts(raw, circ.num_qubits)

    return RunTarget(run_fn=run_fn, noise_model=noise,
                     name=str(name), meta={"stack": "qiskit-ibm-runtime"})


def braket_device(device) -> RunTarget:
    """Wrap an Amazon Braket device (AwsDevice) into a RunTarget via the
    QASM bridge.  Best-effort: exact kwargs vary across Braket SDK
    versions; failures surface as RuntimeError from run_fn."""
    from .noise import NoiseModel
    name = getattr(device, "name", None) or str(device)
    # Braket exposes device calibrations in its own schema; a full mapping
    # lands with a hardware-verified run.  Defaults keep the pipeline usable.
    noise = NoiseModel(getattr(device, "qubit_count", None) or 8)

    def run_fn(circ, *, seed: int = 0, shots: int = 1024) -> dict:
        try:
            from braket.circuits import Circuit as BKCircuit
        except Exception as e:
            raise RuntimeError(
                "amazon-braket-sdk not installed; pip install "
                "amazon-braket-sdk to use braket_device().") from e
        from .io_qasm import to_qasm
        bk = BKCircuit.from_qasm(to_qasm(circ))
        task = device.run(bk, shots=shots)
        result = task.result()
        raw = result.measurement_counts()
        return _reverse_counts(raw, circ.num_qubits)

    return RunTarget(run_fn=run_fn, noise_model=noise,
                     name=str(name), meta={"stack": "amazon-braket-sdk"})
