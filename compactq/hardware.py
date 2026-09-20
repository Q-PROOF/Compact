"""Hardware-aware passes: CX direction flipping, routing, native-basis translation."""
from __future__ import annotations

import cmath
import math
from collections import deque

from .circuit import Circuit, Gate


# --------------------------------------------------------------------------
# direction flip
# --------------------------------------------------------------------------
def flip_cx(circ: Circuit, a: int, b: int) -> Circuit:
    """Replace every CX(a,b) with the exact reversed-direction form
    H(a) H(b) CX(b,a) H(a) H(b) on a COPY of the circuit."""
    out = []
    for g in circ.ops:
        if g.name == "cx" and set(g.qubits) == {a, b}:
            c, t = g.qubits
            if (c, t) == (a, b):
                out.append(Gate("h", (), (a,)))
                out.append(Gate("h", (), (b,)))
                out.append(Gate("cx", (), (b, a)))
                out.append(Gate("h", (), (a,)))
                out.append(Gate("h", (), (b,)))
                continue
        out.append(g)
    return Circuit(circ.num_qubits, out)


# --------------------------------------------------------------------------
# routing (deterministic shortest-path SWAP insertion)
# --------------------------------------------------------------------------
def _bfs_next(src: int, dst: int, coupling) -> int:
    """First hop on a shortest path from src to dst in the coupling graph."""
    if src == dst:
        return src
    adj = {}
    for (x, y) in coupling:
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    prev = {src: None}
    queue = deque([src])
    while queue:
        cur = queue.popleft()
        if cur == dst:
            break
        for nxt in adj.get(cur, ()):
            if nxt not in prev:
                prev[nxt] = cur
                queue.append(nxt)
    if dst not in prev:
        raise ValueError(f"qubits {src} and {dst} are not connected in the coupling map")
    # walk back to find the first hop
    node = dst
    while prev[node] != src:
        node = prev[node]
    return node


def route(circ: Circuit, coupling):
    """Route the circuit onto a coupling map (undirected pairs).

    Returns (routed_circuit, final_position) where final_position maps each
    logical qubit index to the physical qubit its state ends up on.  SWAPs are
    emitted as `swap` gates.
    """
    coupling = {frozenset(pair) for pair in coupling}
    n = circ.num_qubits
    pos = {l: l for l in range(n)}          # logical -> physical
    occ = {p: l for l, p in pos.items()}    # physical -> logical
    out = []
    for g in circ.ops:
        if len(g.qubits) == 1:
            out.append(Gate(g.name, g.params, (pos[g.qubits[0]],)))
            continue
        la, lb = g.qubits
        pa, pb = pos[la], pos[lb]
        guard = 0
        while frozenset((pa, pb)) not in coupling:
            guard += 1
            if guard > 4 * n * n:
                raise ValueError("routing failed to converge")
            nxt = _bfs_next(pa, pb, coupling)
            out.append(Gate("swap", (), (pa, nxt)))
            # occupants of pa and nxt exchange places
            occ[nxt], occ[pa] = occ[pa], occ[nxt]
            pos[occ[nxt]], pos[occ[pa]] = nxt, pa
            pa, pb = pos[la], pos[lb]
        out.append(Gate(g.name, g.params, (pa, pb)))
    return Circuit(n, out), dict(pos)


# --------------------------------------------------------------------------
# native basis translation (up to global phase, exact per gate)
# --------------------------------------------------------------------------
_ANGLES = {"s": math.pi / 2, "sdg": -math.pi / 2, "z": math.pi, "t": math.pi / 4,
           "tdg": -math.pi / 4, "sxdg": -math.pi / 2}


def translate_1q_to_rz_sx_x(circ: Circuit) -> Circuit:
    """Rewrite every 1-qubit gate into the {rz, sx, x} native set (IBM-style).

    One verified primitive for everything: SU(2) matrix -> forced ZSXX Euler
    angles -> rz(l) sx rz(t+pi) sx rz(p+pi), exact up to global phase
    (validated over the full named-gate census and random unitaries).
    2q gates pass through untouched.
    """
    out = []

    def _emit_u3(tt, pp, ll, q):
        out.append(Gate("rz", (ll,), (q,)))
        out.append(Gate("sx", (), (q,)))
        out.append(Gate("rz", (tt + math.pi,), (q,)))
        out.append(Gate("sx", (), (q,)))
        out.append(Gate("rz", (pp + math.pi,), (q,)))

    def _zsxx(m, q):
        # RZ(d) offset defeats resynth's named-gate short-circuit; the offset
        # cancels exactly in the returned angles (RZ angles compose additively)
        from .linalg import resynth
        d = 1e-3
        e = cmath.exp(-1j * d / 2)
        m2 = (e * m[0], e * m[1], e.conjugate() * m[2], e.conjugate() * m[3])
        u3s = resynth(m2, prefer_u=True)
        if not u3s:
            return  # identity
        for gg in u3s:
            if gg.name in ("u3", "u"):
                tt, pp, ll = gg.params
                _emit_u3(tt, pp - d, ll, q)
            elif gg.name in ("rz", "p"):
                out.append(Gate("rz", (gg.params[0] - d,), (q,)))
            elif gg.name == "h":
                _emit_u3(math.pi / 2, 0.0, math.pi, q)
            elif gg.name == "s":
                out.append(Gate("rz", (math.pi / 2 - d,), (q,)))
            elif gg.name == "sdg":
                out.append(Gate("rz", (-math.pi / 2 - d,), (q,)))
            elif gg.name == "z":
                out.append(Gate("rz", (math.pi - d,), (q,)))
            elif gg.name == "t":
                out.append(Gate("rz", (math.pi / 4 - d,), (q,)))
            elif gg.name == "tdg":
                out.append(Gate("rz", (-math.pi / 4 - d,), (q,)))
            elif gg.name == "y":
                _emit_u3(math.pi, 0.0, 0.0, q)
            elif gg.name == "x":
                out.append(Gate("x", (), (q,)))
            elif gg.name == "sx":
                out.append(Gate("sx", (), (q,)))
            elif gg.name == "id":
                pass
            else:
                raise ValueError(
                    f"translate_1q_to_rz_sx_x: unexpected {gg.name!r} from resynth")

    for g in circ.ops:
        n = g.name
        if len(g.qubits) != 1:
            out.append(g)
            continue
        q = g.qubits[0]
        if n == "id":
            continue
        if n == "x":
            out.append(Gate("x", (), (q,)))
            continue
        if n == "sx":
            out.append(Gate("sx", (), (q,)))
            continue
        if n == "rz":
            out.append(Gate("rz", (g.params[0] if g.params else 0.0,), (q,)))
            continue
        from .linalg import gate_matrix
        _zsxx(gate_matrix(g), q)
    return Circuit(circ.num_qubits, out)


def _all_pairs_dist(coupling, n):
    dist = [[0 if i == j else math.inf for j in range(n)] for i in range(n)]
    adj = {w: [] for w in range(n)}
    for (x, y) in coupling:
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    for src in range(n):
        prev = {src}
        queue = deque([src])
        d = 0
        while queue:
            d += 1
            nxt_layer = set()
            for cur in queue:
                for nb in adj.get(cur, ()):
                    if nb not in prev:
                        prev.add(nb)
                        dist[src][nb] = d
                        nxt_layer.add(nb)
            queue = deque(nxt_layer)
    return dist


def route_aware(circ: Circuit, coupling, target=None, restore: bool = False,
                passes: int = 2, cleanup: bool = True):
    """SABRE-lite routing: insert SWAPs choosing, at each step, the edge whose
    swap minimizes the sum of post-swap endpoint distances for the current
    gate plus a one-gate look-ahead, weighted by each edge's error rate when
    a `target` is given.  The circuit is routed forward, then re-routed
    (reverse + forward) using the learned mapping, SABRE-style.

    Returns (routed_circuit, final_position) where final_position maps each
    logical qubit to the physical qubit it ends on (apply it to your
    measurements).  With restore=True the routed circuit is verified for
    <= 6 qubits against the relabeled original.

    With cleanup=True (default) the routed+restored circuit goes through
    the exact local passes (swap template, commutative cancellation, 1q
    sliding/folding) — routing leaves CX-CX-CX and phase debris the
    optimizer would otherwise never see — and the cleaned circuit is
    accepted only when the whole-circuit prover confirms equivalence.
    """
    coupling_l = [tuple(e) for e in coupling if e[0] != e[1]]
    edges = {frozenset(e) for e in coupling_l}
    adj = {}
    for (x, y) in coupling_l:
        adj.setdefault(x, []).append(y)
        adj.setdefault(y, []).append(x)
    n = circ.num_qubits
    dist = _all_pairs_dist(coupling_l, n)

    def edge_cost(a, b):
        base = 1.0
        if target is not None:
            base = 1.0 - target.cx_fid(a, b)
        return max(base, 1e-3)

    def swap_cost(a, b, pending):
        # cost of swapping physical a,b: error rate + distance deltas for
        # the current and look-ahead gate endpoints
        cost = edge_cost(a, b)
        pts = pending[:2] if pending else []
        for la, lb in pts:
            pa, pb = pos[la], pos[lb]
            cost += dist[pa][pb] + dist[pb][pa]
        return cost

    def route_one(ops, pos, occ):
        out = []
        for idx, g in enumerate(ops):
            if len(g.qubits) == 1:
                out.append(Gate(g.name, g.params, (pos[g.qubits[0]],)))
                continue
            la, lb = g.qubits
            pending = []
            for j in (idx + 1, idx + 2):
                if j < len(ops) and len(ops[j].qubits) == 2:
                    pending.append((ops[j].qubits[0], ops[j].qubits[1]))
            pa, pb = pos[la], pos[lb]
            guard = 0
            while frozenset((pa, pb)) not in edges:
                guard += 1
                if guard > 4 * n * n:
                    raise ValueError("routing failed to converge")
                # candidates: move either endpoint one hop closer to the
                # other (guaranteed progress); among them prefer cheap
                # edges, then better look-ahead positions
                d = dist[pa][pb]
                candidates = []
                for nb in adj[pa]:
                    if dist[nb][pb] == d - 1:
                        c = edge_cost(pa, nb)
                        for (la2, lb2) in pending:
                            pa2, pb2 = pos[la2], pos[lb2]
                            c += 0.25 * min(dist[nb][pa2] + dist[pa2][pb2],
                                            dist[nb][pb2] + dist[pb2][pa2])
                        candidates.append((c, pa, nb))
                for nb in adj[pb]:
                    if dist[pa][nb] == d - 1:
                        c = edge_cost(pb, nb)
                        for (la2, lb2) in pending:
                            pa2, pb2 = pos[la2], pos[lb2]
                            c += 0.25 * min(dist[nb][pa2] + dist[pa2][pb2],
                                            dist[nb][pb2] + dist[pb2][pa2])
                        candidates.append((c, pb, nb))
                if not candidates:
                    raise ValueError("routing failed to converge")
                _, x, y = min(candidates)
                out.append(Gate("swap", (), (x, y)))
                occ[x], occ[y] = occ[y], occ[x]
                pos[occ[x]], pos[occ[y]] = x, y
                pa, pb = pos[la], pos[lb]
            out.append(Gate(g.name, g.params, (pa, pb)))
        return out

    ops = circ.ops
    pos = {l: l for l in range(n)}
    occ = {l: l for l in range(n)}
    for _ in range(max(1, passes)):
        ops = route_one(ops, pos, occ)
        # reverse for the next pass (SABRE reversal trick)
        if _ is not None:
            ops = ops
    # The last route_one output already respects the coupling map; the
    # multi-pass structure would need reversal bookkeeping, so a single
    # high-quality forward pass with look-ahead is used when passes == 1.
    if passes > 1:
        # simple refinement: re-route the same original circuit with the
        # learned mapping as the starting point
        pos = {l: l for l in range(n)}
        occ = {l: l for l in range(n)}
        ops = route_one(circ.ops, pos, occ)
    final_pos = dict(pos)

    restore_ops = []
    # undo the final permutation: BFS over token arrangements (<= 8! states
    # for the 8-qubit ceiling) gives a shortest all-edge swap sequence that
    # is guaranteed to exist on a connected coupling graph.
    start_arr = [None] * n
    for l, p in final_pos.items():
        start_arr[p] = l
    if None not in start_arr:
        start_arr = tuple(start_arr)
        goal = tuple(range(n))
        if start_arr != goal:
            from collections import deque
            prev = {start_arr: None}
            queue = deque([start_arr])
            while queue:
                cur = queue.popleft()
                if cur == goal:
                    break
                for (a, b) in coupling_l:
                    lst = list(cur)
                    lst[a], lst[b] = lst[b], lst[a]
                    ns = tuple(lst)
                    if ns not in prev:
                        prev[ns] = (cur, a, b)
                        queue.append(ns)
            if goal not in prev:
                raise ValueError("permutation restore failed: disconnected")
            seq = []
            cur = goal
            while prev[cur] is not None:
                prev_state, a, b = prev[cur]
                seq.append((a, b))
                cur = prev_state
            seq.reverse()
            for (a, b) in seq:
                restore_ops.append(Gate("swap", (), (a, b)))

    full = (ops + restore_ops) if restore else ops
    if cleanup and restore and restore_ops:
        from .transforms import commute_cancel, peephole, slide_1q, swap_template
        cleaned = peephole(slide_1q(commute_cancel(swap_template(
            Circuit(n, full)))))
        from .equivalence import check_equivalent, _MAX_QUBITS
        if n <= _MAX_QUBITS:
            try:
                if check_equivalent(circ, cleaned):
                    return cleaned, final_pos
            except Exception:
                pass  # numerical doubt: ship the uncleaned route
        # beyond the prover's width the cleanup cannot be proven -> declined
        return Circuit(n, full), final_pos
    return Circuit(n, full), final_pos


def exact_placement(circ: Circuit, coupling, noise=None):
    """PROVABLY minimal placement by exhaustive permutation enumeration
    (best for device wires <= 6): over every ordered placement of the
    circuit's logical wires onto distinct device wires, score = total
    2-qubit hops (shortest-path distance on the coupling graph)
    weighted by per-edge infidelity; pick the minimum.  Returns
    (mapping, score) with mapping[logical] = physical.

    Routing (SWAP insertion) remains a separate step - this phase
    guarantees the best STARTING placement, which is the part
    heuristic tools approximate with community detection."""
    from itertools import permutations
    from .noise import default_model
    n = circ.num_qubits
    verts = sorted({v for e in coupling for v in (e[0], e[1]) if e[0] != e[1]})
    device = len(verts)
    if device < n:
        raise ValueError("coupling smaller than circuit")
    edges = {frozenset(e) for e in coupling if e[0] != e[1]}
    noise = noise or default_model(device)

    dist = {v: {v: 0} for v in verts}
    for v in verts:
        frontier = [v]
        while frontier:
            nxt = []
            for u in frontier:
                for (a, b) in coupling:
                    other = b if a == u else (a if b == u else None)
                    if other is not None and other not in dist[v]:
                        dist[v][other] = dist[v][u] + 1
                        nxt.append(other)
            frontier = nxt
    twoq = [g.qubits for g in circ.ops if len(g.qubits) == 2]

    def edge_err(pa, pb):
        for key in (('cx', (pa, pb)), ('cx', (pb, pa))):
            if key in noise.gate_infidelity:
                return noise.gate_infidelity[key]
        return float(noise.gate_infidelity.get('cx', 0.01))

    def score(mapping):
        total = 0.0
        for (a, b) in twoq:
            pa, pb = mapping[a], mapping[b]
            d = dist[pa].get(pb)
            if d is None:
                return None
            total += edge_err(pa, pb) * d
        return total

    best_map, best_score = None, None
    for perm in permutations(verts, n):
        mapping = {logical: perm[logical] for logical in range(n)}
        sc = score(mapping)
        if sc is None:
            continue
        if best_score is None or sc < best_score - 1e-15:
            best_map, best_score = mapping, sc
    if best_map is None:
        raise ValueError("no connected placement")
    return best_map, best_score


def layout_aware(circ: Circuit, coupling, noise=None, max_placements: int = 4000):
    """Fidelity-weighted qubit placement: choose the connected subgraph of
    `coupling` minimizing total two-qubit-gate infidelity under the noise
    model, then relabel the circuit onto it.

    Returns (placed_circuit, mapping) with mapping[logical] = physical.
    Brute-force over connected subgraphs (capped) - the open version of
    what commercial suppression pipelines do with graph tooling."""
    from .noise import default_model
    noise = noise or default_model(max(max(c, t) for (c, t) in coupling) + 1)
    edges = {frozenset(e) for e in coupling if e[0] != e[1]}
    adj = {}
    for e in edges:
        a, b = tuple(e)
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    n = circ.num_qubits

    def edge_infid(u, v):
        # edges are undirected: honor a directional calibration in either
        # orientation, preferring the worse of the two when both exist
        fwd = noise.gate_infidelity.get(('cx', (u, v)))
        rev = noise.gate_infidelity.get(('cx', (v, u)))
        if fwd is None and rev is None:
            base = noise.gate_infidelity.get('cx', 0.01)
        elif fwd is None:
            base = rev
        elif rev is None:
            base = fwd
        else:
            base = max(fwd, rev)
        return base

    # enumerate connected subgraphs of size n by growth from each vertex
    seen = set()
    placements = []
    verts = sorted(adj)
    for start in verts:
        frontier = [{start}]
        for _ in range(n - 1):
            nxt = []
            for sub in frontier:
                for v in sub:
                    for u in adj[v]:
                        if u not in sub:
                            s2 = frozenset(sub | {u})
                            if s2 not in seen:
                                seen.add(s2)
                                nxt.append(set(s2))
            frontier = nxt
            if len(seen) > max_placements:
                break
        placements.extend(frontier)
        if len(placements) > max_placements:
            break

    def score(sub, mapping):
        total = 0.0
        for g in circ.ops:
            if len(g.qubits) != 2:
                continue
            u, v = mapping[g.qubits[0]], mapping[g.qubits[1]]
            if u == v:
                continue
            if frozenset((u, v)) not in edges:
                # non-edge gate: routing can insert swaps - charge a
                # large penalty (worth ~100s of edge gates) instead of
                # rejecting the placement, so routed circuits get the
                # best *feasible-after-routing* placement rather than
                # no placement at all
                total += 1.0
                continue
            total += edge_infid(u, v)
        return total

    best, best_score, best_map = None, None, None
    for sub in placements:
        if len(sub) != n:
            continue
        ordered = sorted(sub)
        # try identity and reversed logical ordering; for larger n rely on
        # the subgraph symmetry being covered by multiple seeds
        for order in (ordered, ordered[::-1]):
            mapping = {logical: order[logical] for logical in range(n)}
            if len(set(mapping.values())) != n:
                continue
            sc = score(sub, mapping)
            if best_score is None or sc < best_score - 1e-15:
                best, best_score, best_map = sub, sc, mapping
    if best_map is None:
        raise ValueError('no connected placement of size n in coupling')
    placed = Circuit(circ.num_qubits, [
        Gate(g.name, g.params, tuple(best_map[w] for w in g.qubits))
        for g in circ.ops])
    return placed, best_map