"""Device metrics in the industry's own language - layer fidelity,
EPLG and 2-qubit-layer structure - estimated from a NoiseModel and
reported BEFORE and AFTER suppression, so the improvement shows up in
the same units hardware vendors report.

Definitions follow the literature honestly but conservatively:
  - a "layer" here is one ASAP layer of 2-qubit gates of the circuit
    (the circuit's own layering, not the device's coupling pattern);
  - the fidelity of one 2-qubit layer is the product of (1 - error)
    over its edges, using the model's per-pair calibration with the
    width fallback of NoiseModel.gate_error;
  - layer fidelity (LF) over L layers is the product of the L layer
    fidelities; EPLG (error per layered gate) is 1 - LF^(1/L).
These are predictions from the calibration data - the same numbers a
vendor would quote for the raw device, plus the suppressed circuit's
 prediction, which no closed service reports.
"""
from __future__ import annotations

from .circuit import Circuit
from .noise import NoiseModel


def two_q_layers(circ: Circuit) -> list[list[tuple[int, int]]]:
    """ASAP layering of the circuit's 2-qubit gates: two gates share a
    layer only when their wire sets are disjoint (1q gates are not
    part of a 2q layer)."""
    occupied_until = {}   # wire -> layer index that occupies it
    layers: list[list[tuple[int, int]]] = []
    for g in circ.ops:
        if len(g.qubits) != 2:
            continue
        a, b = g.qubits
        start = max(occupied_until.get(a, 0), occupied_until.get(b, 0))
        if start == len(layers):
            layers.append([])
        layers[start].append((a, b))
        occupied_until[a] = start + 1
        occupied_until[b] = start + 1
    return layers


def _edge_error(a: int, b: int, noise: NoiseModel) -> float:
    """CX infidelity of edge (a, b) with undirected fallback."""
    key = ('cx', (a, b))
    if key in noise.gate_infidelity:
        return noise.gate_infidelity[key]
    rev = ('cx', (b, a))
    if rev in noise.gate_infidelity:
        return noise.gate_infidelity[rev]
    return float(noise.gate_infidelity.get('cx', 0.008))


def twoq_layer_fidelity(noise: NoiseModel, pairs) -> float:
    """Fidelity of one 2-qubit layer: product of (1 - error) over the
    layer's edges (undirected lookup)."""
    fid = 1.0
    for (a, b) in pairs:
        fid *= 1.0 - _edge_error(a, b, noise)
    return fid


def layer_fidelity(noise: NoiseModel, circ: Circuit,
                   max_layers: int | None = None) -> float:
    """Layer fidelity of `circ` under `noise`: product of the 2q-layer
    fidelities over the first `max_layers` layers (default: all)."""
    lf = 1.0
    layers = two_q_layers(circ)
    if max_layers is not None:
        layers = layers[:max_layers]
    for pairs in layers:
        lf *= twoq_layer_fidelity(noise, pairs)
    return lf


def eplg(noise: NoiseModel, circ: Circuit,
         n_layers: int | None = None) -> float:
    """Error Per Layered Gate: 1 - LF^(1/L).  Returns 0.0 when the
    circuit has no 2-qubit layers."""
    layers = two_q_layers(circ)
    if n_layers is None:
        n_layers = len(layers)
    if n_layers <= 0:
        return 0.0
    lf = layer_fidelity(noise, circ, max_layers=n_layers)
    return 1.0 - lf ** (1.0 / n_layers)


def suppression_metrics(noise: NoiseModel, circ_before: Circuit,
                        circ_after: Circuit) -> dict:
    """Predicted device metrics before and after suppression: layer
    fidelities, EPLG, and the gate/depth counts."""
    lb = layer_fidelity(noise, circ_before)
    la = layer_fidelity(noise, circ_after)
    nb, na = len(two_q_layers(circ_before)), len(two_q_layers(circ_after))
    return {
        "layer_fidelity_before": round(lb, 6),
        "layer_fidelity_after": round(la, 6),
        "layer_fidelity_gain": round(la / lb, 4) if lb > 0 else None,
        "eplg_before": round(eplg(noise, circ_before), 6),
        "eplg_after": round(eplg(noise, circ_after), 6),
        "two_q_layers_before": nb,
        "two_q_layers_after": na,
        "gates_before": len(circ_before.ops),
        "gates_after": len(circ_after.ops),
    }
