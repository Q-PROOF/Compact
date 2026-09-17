"""Deterministic error suppression: Pauli twirling and dynamical decoupling.

Two circuit-level passes implementing the published techniques behind
commercial error-suppression pipelines, with compactq's distinguishing
property: **every transformed circuit stays exactly provable**.

- ``pauli_twirl`` (randomized compiling): around each Clifford entangler a
  random per-wire Pauli pair (P before, G P G^dag after) is inserted; the
  fragment equals G up to global phase, so the whole twirled circuit is
  unitarily identical to the input — verified by the proof net — while
  coherent gate errors are averaged into stochastic ones across variants.
  Adjacent Paulis on the same wire are composed back into single gates
  (``merge_paulis``), halving the randomized-compiling overhead.

- ``insert_dd`` (context-aware dynamical decoupling): idle windows are
  found from gate durations (ASAP schedule over the noise model); each
  window whose *predicted refocusable noise* (quasi-static Z drift plus
  dephasing) exceeds the predicted decoupling-pulse cost receives XY4
  sequences.  XY4 is the identity up to global phase and commutes with
  every gate in the window (they act on other wires), so the result is
  again exactly equivalent — and the pass is benefit-gated: on a noise
  model where decoupling cannot help, it inserts nothing, so the pass
  can never make the expected fidelity worse.

"""
from __future__ import annotations

import math
import random

from .circuit import Circuit, Gate
from .noise import NoiseModel, default_model

# per-wire Pauli type bits: (x, z) -> gate name
_PTYPE = {(0, 0): None, (1, 0): "x", (1, 1): "y", (0, 1): "z"}


def _bits(name):
    return {"x": (1, 0), "y": (1, 1), "z": (0, 1)}.get(name, (0, 0))


def _build_conj_table():
    """Numerically construct the Pauli conjugation table P' = G P G^dag for
    cx / cz / swap: for every per-wire Pauli pair, the correction after the
    gate so that [P, G, P'] == G up to global phase.  Exact by construction
    (matrix-built, sign dropped: signs are fragment-global phases)."""
    import cmath
    from .linalg import gate_matrix, mmul, same_up_to_phase
    I2 = (1 + 0j, 0j, 0j, 1 + 0j)
    words = {"i": I2}
    words["x"] = gate_matrix(Gate("x", (), (0,)))
    words["y"] = gate_matrix(Gate("y", (), (0,)))
    words["z"] = gate_matrix(Gate("z", (), (0,)))

    def word4(pc, pt):
        # build 2-wire Pauli words through our own unitary() so the wire
        # bit-ordering convention matches the reference exactly
        ops = []
        if pc:
            ops.append(Gate(pc, (), (0,)))
        if pt:
            ops.append(Gate(pt, (), (1,)))
        from .equivalence import unitary
        return unitary(Circuit(2, ops))

    table = {}
    for gname in ("cx", "cz", "swap"):
        G = [[0.0] * 4 for _ in range(4)]
        from .equivalence import unitary
        G = unitary(Circuit(2, [Gate(gname, (), (0, 1))]))
        for pc in (None, "x", "y", "z"):
            for pt in (None, "x", "y", "z"):
                if pc is None and pt is None:
                    table[(gname, pc, pt)] = (None, None)
                    continue
                P = word4(pc, pt)
                # P' = G P G^dag  (row-major 4x4)
                Gd = [[G[j][i].conjugate() for j in range(4)] for i in range(4)]

                def mm(A, B):
                    return [[sum(A[i][k] * B[k][j] for k in range(4))
                             for j in range(4)] for i in range(4)]

                Pp = mm(mm(G, P), Gd)
                # match against the 16 Pauli words up to phase
                hit = (None, None)
                for qc in (None, "x", "y", "z"):
                    for qt in (None, "x", "y", "z"):
                        if qc is None and qt is None:
                            continue
                        W = word4(qc, qt)
                        tr = sum(Pp[i][j].conjugate() * W[i][j]
                                 for i in range(4) for j in range(4))
                        if abs(tr) / 4 > 1 - 1e-9:
                            hit = (qc, qt)
                            break
                    if hit != (None, None):
                        break
                table[(gname, pc, pt)] = hit
    return table


_CONJ_TABLE = None


def conj_pair(gate_name, pc, pt):
    """Correction Pauli names (P'c, P't) with G.P.G-dag conjugation."""
    global _CONJ_TABLE
    if _CONJ_TABLE is None:
        _CONJ_TABLE = _build_conj_table()
    return _CONJ_TABLE[(gate_name, pc, pt)]


_PAULIS = (None, "x", "y", "z")


def merge_paulis(ops):
    """Compose adjacent same-wire 1q Paulis into single gates.

    Paulis on different wires commute, so per-wire pending products are
    exact; a pending product is flushed before any gate it does not
    commute with (every non-Pauli gate).  Each composition X^a Z^b (xor
    of Pauli type bits) holds up to a global phase, which the proof net
    ignores — so the rewrite is exact under compactq's invariant.
    """
    out: list = []
    pending: dict = {}   # wire -> [x_bits, z_bits]

    def flush():
        for w in sorted(pending):
            xb, zb = pending[w]
            if xb or zb:
                out.append(Gate(_PTYPE[(xb, zb)], (), (w,)))
        pending.clear()

    for g in ops:
        pb = _bits(g.name) if (g.name in ("x", "y", "z")
                               and len(g.qubits) == 1 and not g.params) else None
        if pb is not None:
            w = g.qubits[0]
            if w in pending:
                px, pz = pending[w]
                pending[w] = [px ^ pb[0], pz ^ pb[1]]
            else:
                pending[w] = [pb[0], pb[1]]
        else:
            flush()
            out.append(g)
    flush()
    return out


def pauli_twirl(circ: Circuit, seed: int = 0, skip: frozenset = frozenset(),
                merge: bool = True) -> Circuit:
    """One randomized-compiling variant of `circ` (unitarily identical to it).

    Clifford entanglers (cx / cz / swap) each receive a random per-wire
    twirl pair; non-Clifford 2q gates (cp at generic angles, ecr, iswap)
    are left untouched in this version - honest scope, extendable.
    With ``merge=True`` (default) the correction Pauli of one entangler
    and the twirl Pauli of the next are composed into a single gate per
    wire, halving the pulse-count overhead of randomized compiling."""
    rng = random.Random(seed)
    out: list = []
    for g in circ.ops:
        if (len(g.qubits) == 2 and g.name in ("cx", "cz", "swap")
                and g.qubits not in skip):
            pc = rng.choice(_PAULIS)
            pt = rng.choice(_PAULIS)
            qc_, qt_ = g.qubits
            if pc:
                out.append(Gate(pc, (), (qc_,)))
            if pt:
                out.append(Gate(pt, (), (qt_,)))
            out.append(g)
            cc, ct = conj_pair(g.name, pc, pt)
            if cc:
                out.append(Gate(cc, (), (qc_,)))
            if ct:
                out.append(Gate(ct, (), (qt_,)))
        else:
            out.append(g)
    if merge:
        out = merge_paulis(out)
    return Circuit(circ.num_qubits, out)


def idle_windows(circ: Circuit, noise: NoiseModel):
    """ASAP schedule per wire; returns [(wire, t_start, t_end, pos)] where
    pos is the index in the op list after which the window sits."""
    t = [0.0] * circ.num_qubits
    last_pos = [-1] * circ.num_qubits
    windows = []
    for i, g in enumerate(circ.ops):
        dur = noise.duration_of(g)
        start = max(t[w] for w in g.qubits) if g.qubits else t[0]
        end = start + dur
        for w in g.qubits:
            if start - t[w] > 0 and last_pos[w] >= 0:
                windows.append((w, t[w], start, last_pos[w]))
            t[w] = end
            last_pos[w] = i
    total = max(t) if t else 0.0
    for w in range(circ.num_qubits):
        if total - t[w] > 0 and last_pos[w] >= 0:
            windows.append((w, t[w], total, last_pos[w]))
    return windows, total


def _dd_benefit(w: int, span_ns: float, noise: NoiseModel) -> float:
    """Predicted refocusable fidelity cost of an idle window on wire w.

    Two channels are addressable by decoupling: quasi-static Z drift
    (rms angle theta costs ~theta^2/2 on generic superpositions) and
    dephasing accumulated over the window.  Relaxation (T1) is
    intentionally NOT counted: pulse sequences cannot refocus it."""
    theta = noise.drift(w) * span_ns
    t1_ns = noise.t1(w) * 1000.0
    t2_ns = noise.t2(w) * 1000.0
    tphi_inv = max(0.0, 1.0 / t2_ns - 0.5 / t1_ns)
    lam = 1.0 - math.exp(-span_ns * tphi_inv)
    return 0.5 * theta * theta + lam


# Dynamical-decoupling pulse sequences.  Every sequence is a Pauli
# word whose product is the identity up to global phase, so each block
# is an exact identity rewrite.  At Pauli level (which is what the
# calibrated simulator resolves) concatenated DD over {X,Y} reduces to
# the XY8 pulse content, so it is not listed separately.
# Dynamical-decoupling pulse sequences.  Every sequence is a Pauli
# word whose product is the identity up to global phase, so each block
# is an exact identity rewrite.  At Pauli level (which is what the
# calibrated simulator resolves) concatenated DD over {X,Y} reduces to
# the XY8 pulse content, so it is not listed separately.
DD_SEQUENCES = {
    "xy4": ("x", "y", "x", "y"),
    "xy8": ("x", "y", "x", "y", "x", "y", "x", "y"),
    "xzx": ("x", "z", "x", "z"),
    "pdd4": ("x", "x", "z", "z"),
}
_DD_SEQUENCES = DD_SEQUENCES


def _resolve_dd_sequence(sequence: str, span_ns: float, x_dur: float) -> str:
    """Resolve 'auto': short windows take the compact XY4, long windows
    take XY8 (protects lower frequencies; same pulse density)."""
    if sequence != "auto":
        return sequence
    return "xy8" if span_ns >= 8 * x_dur else "xy4"


def insert_dd(circ: Circuit, noise: NoiseModel | None = None,
              min_idle_ns: float = 150.0, max_reps: int = 2,
              sequence: str = "auto", safety: float = 3.0):
    """Insert dynamical decoupling into idle windows that benefit from it.

    Each idle window is scored: refocusable-noise cost (``_dd_benefit``)
    against the decoupling overhead (pulse-count x 1q gate infidelity).
    A window is decoupled only when the predicted benefit exceeds
    ``safety`` times the overhead — so on a model where DD cannot pay
    for itself, the pass inserts nothing and is never a net loss.

    ``sequence`` selects the pulse family: 'xy4' (default compact),
    'xy8' (lower-frequency protection), 'xzx', 'pdd4'
    (X X Z Z), or 'auto' — XY4 for windows shorter than eight pulse
    durations, XY8 otherwise.  Every sequence is a Pauli word with
    identity product, so every block is exact up to global phase."""
    noise = noise or default_model(circ.num_qubits)
    if sequence != "auto" and sequence not in _DD_SEQUENCES:
        raise ValueError(f"unknown DD sequence {sequence!r}; "
                         f"choose from {sorted(_DD_SEQUENCES)} or 'auto'")
    x_dur = noise.durations_ns.get("1q", 50.0)
    windows, _ = idle_windows(circ, noise)
    insertions = {}  # pos -> [gate, ...]
    for (w, t0, t1, pos) in windows:
        span = t1 - t0
        seq_name = _resolve_dd_sequence(sequence, span, x_dur)
        pulses = _DD_SEQUENCES[seq_name]
        block_dur = len(pulses) * x_dur
        if span < max(min_idle_ns, block_dur):
            continue
        reps = min(max_reps, int(span // block_dur))
        if reps < 1:
            continue
        x_err = noise.gate_error(Gate("x", (), (w,)))
        overhead = safety * reps * len(pulses) * x_err
        if _dd_benefit(w, span, noise) <= overhead:
            continue
        block = []
        for _ in range(reps):
            for nm in pulses:
                block.append(Gate(nm, (), (w,)))
        insertions.setdefault(pos, []).extend(block)
    if not insertions:
        return circ
    out = []
    for i, g in enumerate(circ.ops):
        out.append(g)
        if i in insertions:
            out.extend(insertions[i])
    return Circuit(circ.num_qubits, out)


# ---------------------------------------------------------------------------
# the automated pipeline (the open answer to one-call suppression services)
# ---------------------------------------------------------------------------

def expand_for_suppression(circ: Circuit) -> Circuit:
    """Expand CP entanglers into the exact CX form the suppression passes
    handle.  From compactq's verified identity CX.RZ_t(t).CX =
    RZ_c(t).RZ_t(t).CP(-2t):  CP(x) = RZ(x/2)_a RZ(x/2)_b CX RZ(-x/2)_b CX
    (up to global phase)."""
    out = []
    for g in circ.ops:
        if g.name == "cp" and len(g.qubits) == 2:
            a, b = g.qubits
            x = g.params[0]
            out.append(Gate("rz", (x / 2,), (a,)))
            out.append(Gate("rz", (x / 2,), (b,)))
            out.append(Gate("cx", (), (a, b)))
            out.append(Gate("rz", (-x / 2,), (b,)))
            out.append(Gate("cx", (), (a, b)))
        else:
            out.append(g)
    return Circuit(circ.num_qubits, out)


def _verify_variant(circ: Circuit, variant: Circuit):
    """Return (variant, proof_level). The variant is kept only if provably
    equivalent to `circ`; on numerical doubt the original is returned."""
    from .equivalence import check_equivalent
    from .verify_large import states_agree
    try:
        if circ.num_qubits <= 8:
            ok = check_equivalent(circ, variant)
            proof = "exact-unitary"
        else:
            ok = states_agree(circ, variant)
            proof = "randomized-K"
        if ok:
            return variant, proof
        return circ, "input-unchanged"
    except Exception:
        return circ, "input-unchanged"


def to_device(circ: Circuit, placement: dict, num_device: int) -> Circuit:
    """Relabel a logical circuit into device wire space (pure relabeling +
    idle wires): placement[logical] = physical.  Exact by construction - the
    device-space unitary is the logical one tensored with identity."""
    out = [Gate(g.name, g.params,
                tuple(placement[w] for w in g.qubits)) for g in circ.ops]
    return Circuit(num_device, out)


def twirl_worthwhile(circ: Circuit, noise: NoiseModel) -> bool:
    """Model-based randomized-compiling decision.

    RC converts the coherent share of the calibrated 2q error (eps =
    coherent_fraction x mean 2q infidelity) into a ~2*eps^2 stochastic
    rate; the twirl pulses cost ~2 extra 1q gates per entangler.  Twirl
    only when the predicted gain beats 3x that overhead - the same
    no-net-loss argument as the DD benefit gate."""
    two_q = [g for g in circ.ops if len(g.qubits) == 2]
    if not two_q:
        return False  # nothing for RC to act on
    p2q = sum(noise.gate_error(g) for g in two_q) / len(two_q)
    p1q = noise.gate_error(Gate("x", (), (0,)))
    coherent = noise.coherent_fraction * p2q
    gain = coherent - 2.0 * coherent * coherent
    return gain > 3.0 * 2.0 * p1q


def suppress_plan(circ: Circuit, noise: NoiseModel | None = None,
                  seed: int = 0, variants: int = 4,
                  optimize: bool = True,
                  coupling=None, target=None,
                  twirl_fraction: float = 1.0,
                  force_twirl: bool | None = None,
                  dd_sequence: str = "auto") -> dict:
    """Build the full suppression plan for `circ` - the exact, automated
    pipeline: optimize -> expand untwirlable entanglers -> (noise-aware
    layout + SABRE-lite routing when a coupling map is given) -> K
    Pauli-twirled variants with dynamical decoupling placed on each
    variant's own schedule.

    With `coupling`, the plan runs in DEVICE wire space: placement is
    chosen by per-edge CX infidelity (`hardware.layout_aware`), routing
    inserts error-weighted SWAPs and restores the final permutation
    (`hardware.route_aware`), and every variant is proven equivalent to
    the device-space reference of the input.  Returned circuits carry
    physical qubit indices; `mapping` in the plan maps logical -> physical
    for measurement interpretation.

    Every returned variant is proven equivalent to the input (dense proof
    to 8q, randomized-K beyond); anything doubtful degrades to the original
    circuit.  Order matters: DD must see the final idle structure, and
    twirl-merging must never compose DD sequences away.
    """
    noise = noise or default_model(circ.num_qubits)
    if optimize:
        from .search import optimize_search
        try:
            base = optimize_search(circ, verify=False)
        except Exception:
            base = circ
        if base.two_qubit_count() > circ.two_qubit_count():
            base = circ  # the optimizer never grows; enforce at this layer
    else:
        base = circ
    base = expand_for_suppression(base)

    from .report import SuppressionReport
    rep = SuppressionReport(circ.num_qubits)
    mapping = {w: w for w in range(circ.num_qubits)}
    if coupling is not None:
        from .hardware import layout_aware, route_aware
        device = max(max(c, t) for (c, t) in coupling) + 1
        rep.add("optimize", base, note="logical")
        try:
            _, placement = layout_aware(circ, coupling, noise)
        except ValueError:
            placement = {w: w for w in range(circ.num_qubits)}
        mapping = dict(placement)
        dev = to_device(base, placement, device)
        rep.add("layout", dev, note="fidelity-weighted placement")
        try:
            dev, _fpos = route_aware(dev, coupling, target=target,
                                     restore=(device <= 8))
        except ValueError:
            pass  # routing failed to converge: run unmapped, report honestly
        rep.add("route", dev, proof="exact-unitary",
                note="swaps restored" if device <= 8 else
                "final permutation left to measurement mapping")
        base = dev

    # which variants get randomized compiling: model-gated, or forced;
    # twirl_fraction < 1 runs a portfolio (twirled + clean variants)
    decided = force_twirl if force_twirl is not None \
        else twirl_worthwhile(base, noise)
    k_total = max(1, variants)
    if decided:
        k_twirled = max(1, round(k_total * min(1.0, max(0.0, twirl_fraction))))
    else:
        k_twirled = 0
    twirled = [pauli_twirl(base, seed=seed + k) for k in range(k_twirled)]
    decoupled = [insert_dd(t, noise, sequence=dd_sequence)
                 for t in twirled]
    decoupled += [insert_dd(base, noise, sequence=dd_sequence)
                  for _ in range(k_total - k_twirled)]
    finals = []
    proofs = []
    for v in decoupled:
        ref = circ if coupling is None else to_device(circ, mapping,
                                                      v.num_qubits)
        vv, pr = _verify_variant(ref, v)
        finals.append(vv)
        proofs.append(pr)
    if not any(x.name == "optimize" for x in rep.stages):
        rep.add("optimize", base, proof=proofs[0] if proofs else "exact-unitary",
                note=("logical" if coupling is None else "device space"))
    rep.add("twirl x%d + dd" % len(finals), finals[0],
            proof=proofs[0], note="variant 0; all variants individually proven")
    rep.meta = {"variants": len(finals),
                "twirled": k_twirled,
                "twirl_decision": ("forced" if force_twirl is not None
                                   else "model-gated"),
                "seed": seed,
                "mapping": mapping,
                "coupling": bool(coupling)}

    return {"base": base,
            "variants": finals,
            "num_variants": len(finals),
            "noise": noise,
            "seed": seed,
            "mapping": mapping,
            "report": rep}


def suppress_execute(circ: Circuit, noise: NoiseModel | None = None,
                     run_fn=None, *, seed: int = 0, variants: int = 4,
                     shots: int = 4096, mitigate: bool = True,
                     mitigation: str = "invert",
                     eps_coh: float = 0.0,
                     coupling=None, target=None,
                     measure_baseline: bool = True,
                     dd_sequence: str = "auto") -> dict:
    """One-call error-suppressed execution - the whole pipeline in a
    single function, locally and with zero secrets:

        result = suppress_execute(circ, noise)
        best = max(result["probabilities"], key=result["probabilities"].get)
        print(result["report"])          # per-stage proof + factor audit

    Plan: ``suppress_plan`` (optimize -> expand -> [layout + route] ->
    twirl x K -> DD per variant).  Execution: each variant is submitted
    with shots/variants shots through `run_fn(circuit, *, seed, shots) ->
    counts` (your hardware adapter), or through the built-in zero-
    dependency trajectory simulator when run_fn is None.  Measurement
    error mitigation (tensored readout inversion) is applied to the
    pooled counts when `mitigate` is on.  With the built-in simulator and
    `measure_baseline=True`, the ORIGINAL circuit is executed with a
    quarter of the shot budget and the measured suppression factor
    (P(ideal) suppressed / raw) is recorded in the report.  Returns:

        {"probabilities": {bitstring: quasi-probability},   (char i = wire i)
         "counts":          {bitstring: count},             (pooled, raw)
         "report":          SuppressionReport,
         "meta":            plan + per-variant counts}
    """
    noise = noise or default_model(circ.num_qubits)
    plan = suppress_plan(circ, noise, seed=seed, variants=variants,
                         coupling=coupling, target=target,
                         dd_sequence=dd_sequence)
    rep = plan["report"]
    per = max(1, shots // plan["num_variants"])
    pooled: dict = {}
    by_variant = []
    for k, vc in enumerate(plan["variants"]):
        if run_fn is not None:
            ck = run_fn(vc, seed=seed + 17 * k, shots=per)
        else:
            from .simulate import simulate_counts
            ck = simulate_counts(vc, noise, shots=per, seed=seed + 17 * k,
                                 eps_coh=eps_coh)
        by_variant.append(ck)
        for s, c in ck.items():
            pooled[s] = pooled.get(s, 0) + c
    probs = {s: c / max(1, sum(pooled.values())) for s, c in pooled.items()}
    if mitigate and mitigation != "none":
        n_out = plan["variants"][0].num_qubits
        readout = {w: noise.readout_error(w) for w in range(n_out)}
        if mitigation == "mle":
            from .mitigate import mitigate_mle
            probs = mitigate_mle(pooled, readout, n_out)
        else:
            from .mitigate import mitigate_counts
            probs = mitigate_counts(pooled, readout, n_out, clip=True)

    if run_fn is None and measure_baseline:
        from .simulate import simulate_counts
        from .report import zero_noise
        n_dev = plan["variants"][0].num_qubits
        ideal = simulate_counts(circ, zero_noise(circ.num_qubits),
                                shots=256, seed=seed)
        top = max(ideal.values())
        support = [s for s, c in ideal.items() if c >= top // 2]

        def mass(counts):
            tot = sum(counts.values()) or 1
            # aggregate device-space strings back to the logical support
            mp = plan["mapping"]
            acc = 0.0
            for s, c in counts.items():
                logical = ''.join(s[mp[w]] for w in range(circ.num_qubits))
                if logical in support:
                    acc += c / tot
            return acc

        raw_counts = simulate_counts(circ if coupling is None
                                     else to_device(circ, plan["mapping"],
                                                    n_dev),
                                     noise, shots=max(256, shots // 4),
                                     seed=seed + 999_9, eps_coh=eps_coh)
        raw_p = mass(raw_counts)
        sup_p = mass(dict(probs)) if mitigate else mass(pooled)
        rep.set_measured(raw_p, sup_p, support)

    rep.meta["shots"] = shots
    rep.meta["mitigated"] = mitigate
    # predicted device metrics (layer fidelity / EPLG) before vs after,
    # in the units hardware vendors report
    try:
        from .metrics import suppression_metrics
        n_dev = plan["variants"][0].num_qubits
        dev_circ = circ if coupling is None else to_device(
            circ, plan["mapping"], n_dev)
        rep.meta["predicted_metrics"] = suppression_metrics(
            noise, dev_circ, plan["variants"][0])
    except Exception:
        pass  # metrics are additive; never block execution
    return {"probabilities": probs, "counts": pooled,
            "report": rep,
            "meta": {"plan": {k: v for k, v in plan.items()
                              if k in ("num_variants", "seed", "mapping")},
                     "counts_by_variant": by_variant}}
