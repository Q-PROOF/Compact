//! Rust acceleration for compactq's KAK hot path.
//!
//! Zero-dependency core contract: compactq works without this module.  When the
//! compiled extension is present, `compactq.kak` uses it for the two hottest
//! numeric kernels — block-unitary construction and Weyl class decisions —
//! with bit-for-bit-compatible semantics to the Python reference.
//!
//! Matrices are row-major; complex numbers are (re, im) f64 pairs.

pub mod resynth_1q;
use resynth_1q::resynth_1q as py_resynth_1q;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

type Mat = Vec<f64>; // 2*n*n floats: re, im pairs

/// U = A @ B for square complex matrices given as interleaved re/im pairs.
fn matmul(a: &Mat, b: &Mat, n: usize) -> Mat {
    let mut out = vec![0.0f64; 2 * n * n];
    for i in 0..n {
        for k in 0..n {
            let ar = a[2 * (i * n + k)];
            let ai = a[2 * (i * n + k) + 1];
            if ar == 0.0 && ai == 0.0 {
                continue;
            }
            for j in 0..n {
                let br = b[2 * (k * n + j)];
                let bi = b[2 * (k * n + j) + 1];
                out[2 * (i * n + j)] += ar * br - ai * bi;
                out[2 * (i * n + j) + 1] += ar * bi + ai * br;
            }
        }
    }
    out
}

/// Complex Kronecker product a (x) b with a acting on the high bit.
fn kron(a: &Mat, b: &Mat, na: usize, nb: usize) -> Mat {
    let n = na * nb;
    let mut out = vec![0.0f64; 2 * n * n];
    for ia in 0..na {
        for ja in 0..na {
            for ib in 0..nb {
                for jb in 0..nb {
                    let ar = a[2 * (ia * na + ja)];
                    let ai = a[2 * (ia * na + ja) + 1];
                    let br = b[2 * (ib * nb + jb)];
                    let bi = b[2 * (ib * nb + jb) + 1];
                    let i = ia * nb + ib;
                    let j = ja * nb + jb;
                    out[2 * (i * n + j)] = ar * br - ai * bi;
                    out[2 * (i * n + j) + 1] = ar * bi + ai * br;
                }
            }
        }
    }
    out
}

fn single_qubit_matrix(name: &str, p: f64) -> Option<Mat> {
    match name {
        "i" | "id" | "u0" => return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
        "x" => return Some(vec![0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
        "y" => return Some(vec![0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0]),
        "z" => return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0]),
        "h" => {
            let s = std::f64::consts::FRAC_1_SQRT_2;
            return Some(vec![s, 0.0, s, 0.0, s, 0.0, -s, 0.0]);
        }
        "s" => return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        "sdg" => return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0]),
        "t" => {
            let c = std::f64::consts::FRAC_1_SQRT_2;
            return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, c, c]);
        }
        "tdg" => {
            let c = std::f64::consts::FRAC_1_SQRT_2;
            return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, c, -c]);
        }
        "sx" => {
            // e^{i pi/4} RX(pi/2)
            return Some(vec![0.5, 0.5, 0.5, -0.5, 0.5, -0.5, 0.5, 0.5]);
        }
        "sxdg" => {
            return Some(vec![0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5]);
        }
        "rx" => {
            let (c, s) = ((p / 2.0).cos(), (p / 2.0).sin());
            return Some(vec![c, 0.0, 0.0, -s, 0.0, -s, c, 0.0]);
        }
        "ry" => {
            let (c, s) = ((p / 2.0).cos(), (p / 2.0).sin());
            return Some(vec![c, 0.0, -s, 0.0, s, 0.0, c, 0.0]);
        }
        "rz" => {
            let (c, s) = ((p / 2.0).cos(), (p / 2.0).sin());
            return Some(vec![c, -s, 0.0, 0.0, 0.0, 0.0, c, s]);
        }
        "p" => {
            return Some(vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, p.cos(), p.sin()]);
        }
        "u3" | "u" => return None, // 3-param; handled by caller in Python
        _ => return None,
    }
}

/// Build the unitary of a circuit fragment on `n` qubits.
/// ops: flat list of (name, params, qubits).
/// Returns the 2^n x 2^n matrix as interleaved re/im f64 pairs (row-major).
#[pyfunction]
#[pyo3(signature = (n, names, params_flat, param_counts, qubits_flat, qubit_counts))]
fn block_unitary(
    n: usize,
    names: Vec<String>,
    params_flat: Vec<f64>,
    param_counts: Vec<usize>,
    qubits_flat: Vec<usize>,
    qubit_counts: Vec<usize>,
) -> PyResult<Mat> {
    if names.len() != param_counts.len() || names.len() != qubit_counts.len() {
        return Err(PyValueError::new_err("ragged op lists"));
    }
    let dim = 1usize << n;
    let mut u: Mat = vec![0.0; 2 * dim * dim];
    for i in 0..dim {
        u[2 * (i * dim + i)] = 1.0;
    }
    let mut poff = 0usize;
    let mut qoff = 0usize;
    for (oi, name) in names.iter().enumerate() {
        let pc = param_counts[oi];
        let qc = qubit_counts[oi];
        let params = &params_flat[poff..poff + pc];
        let qubits = &qubits_flat[qoff..qoff + qc];
        poff += pc;
        qoff += qc;
        let g: Mat = match (name.as_str(), qc) {
            ("cx", 2) => {
                let (c, t) = (qubits[0], qubits[1]);
                let mut g = vec![0.0; 2 * dim * dim];
                for r in 0..dim {
                    // control set: the row swaps with its target-flipped twin
                    let src = if (r >> c) & 1 == 1 { r ^ (1 << t) } else { r };
                    g[2 * (r * dim + src)] = 1.0;
                }
                g
            }
            ("cz", 2) => {
                let mut g = vec![0.0; 2 * dim * dim];
                for r in 0..dim {
                    let sign = if ((r >> qubits[0]) & 1) == 1 && ((r >> qubits[1]) & 1) == 1 {
                        -1.0
                    } else {
                        1.0
                    };
                    g[2 * (r * dim + r)] = sign;
                }
                g
            }
            ("swap", 2) => {
                let (a, b) = (qubits[0], qubits[1]);
                let mut g = vec![0.0; 2 * dim * dim];
                for r in 0..dim {
                    let ra = (r >> a) & 1;
                    let rb = (r >> b) & 1;
                    let src = if ra != rb { r ^ ((1 << a) | (1 << b)) } else { r };
                    g[2 * (r * dim + src)] = 1.0;
                }
                g
            }
            ("cp", 2) => {
                let ph = if pc == 1 { params[0] } else { 0.0 };
                let mut g = vec![0.0; 2 * dim * dim];
                for r in 0..dim {
                    g[2 * (r * dim + r)] = 1.0;
                }
                let both = (1 << qubits[0]) | (1 << qubits[1]);
                for r in 0..dim {
                    if r & both == both {
                        g[2 * (r * dim + r)] = ph.cos();
                        g[2 * (r * dim + r) + 1] = ph.sin();
                    }
                }
                g
            }
            _ if qc == 1 => {
                let q = qubits[0];
                let m = if name == "u3" || name == "u" {
                    if pc != 3 { return Err(PyValueError::new_err("u3 needs 3 params")); }
                    // RZ(ph) . RY(tt) . RZ(lam) (compactq's u3 semantics)
                    let rz1 = single_qubit_matrix("rz", params[1]).unwrap();
                    let ry = single_qubit_matrix("ry", params[0]).unwrap();
                    let rz2 = single_qubit_matrix("rz", params[2]).unwrap();
                    matmul(&matmul(&rz1, &ry, 2), &rz2, 2)
                } else {
                    match single_qubit_matrix(name, if pc > 0 { params[0] } else { 0.0 }) {
                        Some(m) => m,
                        None => {
                            return Err(PyValueError::new_err(format!(
                                "unsupported gate {name} with {qc} qubits"
                            )))
                        }
                    }
                };
                // embed g on wire q: identity everywhere except the 2x2 block
                let mut gfull = vec![0.0f64; 2 * dim * dim];
                for r in 0..dim {
                    for c in 0..dim {
                        if (r & !(1usize << q)) == (c & !(1usize << q)) {
                            let rq = (r >> q) & 1;
                            let cq = (c >> q) & 1;
                            gfull[2 * (r * dim + c)] = m[2 * (rq * 2 + cq)];
                            gfull[2 * (r * dim + c) + 1] = m[2 * (rq * 2 + cq) + 1];
                        }
                    }
                }
                gfull
            }
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unsupported gate {name} with {qc} qubits"
                )))
            }
        };
        u = matmul(&g, &u, dim); // gate on the left: circuit order
    }
    Ok(u)
}

/// Determinant of a complex 4x4 matrix (interleaved re/im).
#[pyfunction]
fn det4(m: Vec<f64>) -> PyResult<(f64, f64)> {
    if m.len() != 32 {
        return Err(PyValueError::new_err("expected 32 floats (16 complex)"));
    }
    let g = |i: usize, j: usize| (m[2 * (i * 4 + j)], m[2 * (i * 4 + j) + 1]);
    let cmul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let csub = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0, a.1 - b.1);
    let mut det = (0.0, 0.0);
    for j in 0..4 {
        let mut sub = [(0.0f64, 0.0f64); 9];
        let mut idx = 0;
        for r in 1..4 {
            for c in 0..4 {
                if c != j {
                    sub[idx] = g(r, c);
                    idx += 1;
                }
            }
        }
        let t0 = cmul(sub[0], csub(cmul(sub[4], sub[8]), cmul(sub[5], sub[7])));
        let t1 = cmul(sub[1], csub(cmul(sub[3], sub[8]), cmul(sub[5], sub[6])));
        let t2 = cmul(sub[2], csub(cmul(sub[3], sub[7]), cmul(sub[4], sub[6])));
        let d3 = (t0.0 - t1.0 + t2.0, t0.1 - t1.1 + t2.1);
        let term = cmul(g(0, j), d3);
        if j % 2 == 0 {
            det = (det.0 + term.0, det.1 + term.1);
        } else {
            det = (det.0 - term.0, det.1 - term.1);
        }
    }
    Ok(det)
}


/// Build the full 2^n x 2^n unitary of a circuit by SPARSE gate application
/// (row-pair updates for 1q gates, row swaps / phase flips for 2q gates),
/// which stays fast at 7-8 qubits where dense matmul would not.
#[pyfunction]
#[pyo3(signature = (n, names, params_flat, param_counts, qubits_flat, qubit_counts))]
fn sim_unitary(
    n: usize,
    names: Vec<String>,
    params_flat: Vec<f64>,
    param_counts: Vec<usize>,
    qubits_flat: Vec<usize>,
    qubit_counts: Vec<usize>,
) -> PyResult<Mat> {
    if names.len() != param_counts.len() || names.len() != qubit_counts.len() {
        return Err(PyValueError::new_err("ragged op lists"));
    }
    let dim = 1usize << n;
    let mut u: Mat = vec![0.0; 2 * dim * dim];
    for i in 0..dim {
        u[2 * (i * dim + i)] = 1.0;
    }
    let mut poff = 0usize;
    let mut qoff = 0usize;
    for (oi, name) in names.iter().enumerate() {
        let pc = param_counts[oi];
        let qc = qubit_counts[oi];
        let params = &params_flat[poff..poff + pc];
        let qubits = &qubits_flat[qoff..qoff + qc];
        poff += pc;
        qoff += qc;
        match (name.as_str(), qc) {
            ("cx", 2) => {
                let (c, t) = (qubits[0], qubits[1]);
                for r in 0..dim {
                    if (r >> c) & 1 == 1 && (r >> t) & 1 == 0 {
                        let r2 = r | (1 << t);
                        for col in 0..dim {
                            for k in 0..2 {
                                let a = u[2 * (r * dim + col) + k];
                                u[2 * (r * dim + col) + k] = u[2 * (r2 * dim + col) + k];
                                u[2 * (r2 * dim + col) + k] = a;
                            }
                        }
                    }
                }
            }
            ("cz", 2) => {
                let m = (1usize << qubits[0]) | (1usize << qubits[1]);
                for r in 0..dim {
                    if r & m == m {
                        for col in 0..dim {
                            u[2 * (r * dim + col)] = -u[2 * (r * dim + col)];
                            u[2 * (r * dim + col) + 1] = -u[2 * (r * dim + col) + 1];
                        }
                    }
                }
            }
            ("swap", 2) => {
                let (a, b) = (qubits[0], qubits[1]);
                for r in 0..dim {
                    let ra = (r >> a) & 1;
                    let rb = (r >> b) & 1;
                    if ra == 1 && rb == 0 {
                        let r2 = r ^ ((1 << a) | (1 << b));
                        for col in 0..dim {
                            for k in 0..2 {
                                let x = u[2 * (r * dim + col) + k];
                                u[2 * (r * dim + col) + k] = u[2 * (r2 * dim + col) + k];
                                u[2 * (r2 * dim + col) + k] = x;
                            }
                        }
                    }
                }
            }
            ("cp", 2) => {
                let ph = if pc == 1 { params[0] } else { 0.0 };
                let (cr, ci) = (ph.cos(), ph.sin());
                let m = (1usize << qubits[0]) | (1usize << qubits[1]);
                for r in 0..dim {
                    if r & m == m {
                        for col in 0..dim {
                            let re = u[2 * (r * dim + col)];
                            let im = u[2 * (r * dim + col) + 1];
                            u[2 * (r * dim + col)] = re * cr - im * ci;
                            u[2 * (r * dim + col) + 1] = re * ci + im * cr;
                        }
                    }
                }
            }
            _ if qc == 1 => {
                let q = qubits[0];
                let m2 = if name == "u3" || name == "u" {
                    if pc != 3 {
                        return Err(PyValueError::new_err("u3 needs 3 params"));
                    }
                    let rz1 = single_qubit_matrix("rz", params[1]).unwrap();
                    let ry = single_qubit_matrix("ry", params[0]).unwrap();
                    let rz2 = single_qubit_matrix("rz", params[2]).unwrap();
                    matmul(&matmul(&rz1, &ry, 2), &rz2, 2)
                } else {
                    match single_qubit_matrix(name, if pc > 0 { params[0] } else { 0.0 }) {
                        Some(m) => m,
                        None => {
                            return Err(PyValueError::new_err(format!(
                                "unsupported gate {name} with {qc} qubits"
                            )))
                        }
                    }
                };
                let bit = 1usize << q;
                // (m2 @ U) applied to row pairs: rows r (bit=0) and r|bit
                for r in 0..dim {
                    if r & bit != 0 {
                        continue;
                    }
                    let r1 = r | bit;
                    for col in 0..dim {
                        let x0r = u[2 * (r * dim + col)];
                        let x0i = u[2 * (r * dim + col) + 1];
                        let x1r = u[2 * (r1 * dim + col)];
                        let x1i = u[2 * (r1 * dim + col) + 1];
                        let m00r = m2[0]; let m00i = m2[1];
                        let m01r = m2[2]; let m01i = m2[3];
                        let m10r = m2[4]; let m10i = m2[5];
                        let m11r = m2[6]; let m11i = m2[7];
                        u[2 * (r * dim + col)] =
                            m00r * x0r - m00i * x0i + m01r * x1r - m01i * x1i;
                        u[2 * (r * dim + col) + 1] =
                            m00r * x0i + m00i * x0r + m01r * x1i + m01i * x1r;
                        u[2 * (r1 * dim + col)] =
                            m10r * x0r - m10i * x0i + m11r * x1r - m11i * x1i;
                        u[2 * (r1 * dim + col) + 1] =
                            m10r * x0i + m10i * x0r + m11r * x1i + m11i * x1r;
                    }
                }
            }
            _ => {
                return Err(PyValueError::new_err(format!(
                    "unsupported gate {name} with {qc} qubits"
                )))
            }
        }
    }
    Ok(u)
}

/// |Tr(A^dag B)| / dim over two interleaved re/im unitaries of equal dim.
#[pyfunction]
fn trace2(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    if a.len() != b.len() || a.len() % 2 != 0 || a.is_empty() {
        return Err(PyValueError::new_err(
            "trace2 expects two equal-length interleaved re/im buffers",
        ));
    }
    let dim2 = a.len() / 2;
    let dim = (dim2 as f64).sqrt().round() as usize;
    if dim * dim != dim2 {
        return Err(PyValueError::new_err("buffer is not a square matrix"));
    }
    let mut re = 0.0f64;
    let mut im = 0.0f64;
    for i in 0..dim {
        for j in 0..dim {
            let k = 2 * (i * dim + j);
            let ar = a[k];
            let ai = a[k + 1];
            let br = b[k];
            let bi = b[k + 1];
            re += ar * br + ai * bi;
            im += ai * br - ar * bi;
        }
    }
    Ok((re * re + im * im).sqrt() / dim as f64)
}
/// BFS parity-network synthesis (port of compactq.parity._parity_network).
/// Each wire carries an n-bit parity mask; the state packs all n masks into
/// an n*n-bit integer (wire w occupies bits w*n .. w*n+n-1).  terms:
/// interleaved (mask, angle).  Returns (cx_pairs, rz_list) with the same
/// replay convention as the Python reference.
#[pyfunction]
#[pyo3(signature = (n, terms, max_states))]
fn parity_network(n: usize, terms: Vec<(u64, f64)>, max_states: usize) -> PyResult<Vec<(u8, usize, f64)>> {
    if n > 8 {
        return Err(PyValueError::new_err("n > 8 unsupported"));
    }
    let supports: Vec<(u64, f64)> = terms.into_iter().filter(|(s, a)| *s != 0 && a.abs() > 1e-12).collect();
    if supports.is_empty() {
        return Err(PyValueError::new_err("no supports"));
    }
    // state packing: wire w mask at bits w*n..
    let pack = |masks: &[u64]| -> u64 {
        let mut st = 0u64;
        for (w, m) in masks.iter().enumerate() {
            st |= (m & ((1u64 << n) - 1)) << (w * n);
        }
        st
    };
    let wire_mask = |st: u64, w: usize| -> u64 { (st >> (w * n)) & ((1u64 << n) - 1) };
    let start_masks: Vec<u64> = (0..n).map(|w| 1u64 << w).collect();
    let start = pack(&start_masks);
    let covered = |st: u64, sup: &[(u64, f64)]| sup.iter().all(|(s, _)| {
        (0..n).any(|w| wire_mask(st, w) == *s)
    });
        if covered(start, &supports) {
        return Ok(supports.iter().map(|(s, a)| (0u8, (*s).trailing_zeros() as usize, *a)).collect());
    }
    let mut seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    seen.insert(start);
    let mut parents: std::collections::HashMap<u64, (u64, usize, usize)> =
        std::collections::HashMap::new();
    let mut frontier: Vec<u64> = vec![start];
    let mut goal: Option<u64> = None;
    'outer: while !frontier.is_empty() {
        let mut next_frontier: Vec<u64> = vec![];
        for st in frontier.iter() {
            if seen.len() >= max_states { break 'outer; }
            for c in 0..n {
                for t in 0..n {
                    if c == t { continue; }
                    let mc = wire_mask(*st, c);
                    let mt = wire_mask(*st, t);
                    let nmasks: Vec<u64> = (0..n).map(|w| if w == t { mt ^ mc } else { wire_mask(*st, w) }).collect();
                    let ns = pack(&nmasks);
                    if seen.contains(&ns) { continue; }
                    seen.insert(ns);
                    parents.insert(ns, (*st, c, t));
                    if covered(ns, &supports) {
                        goal = Some(ns);
                        break 'outer;
                    }
                    next_frontier.push(ns);
                }
            }
        }
        frontier = next_frontier;
    }
    let goal_state = match goal {
        Some(g) => g,
        None => return Err(PyValueError::new_err("budget exhausted")),
    };
    let mut path_rev: Vec<(usize, usize)> = vec![];
    let mut cur = goal_state;
    while cur != start {
        let (prev, c, t) = parents[&cur];
        path_rev.push((c, t));
        cur = prev;
    }
    path_rev.reverse();
    let path = path_rev;    // replay with per-wire masks (Python-reference convention): one ordered
    // op sequence, kind 0 = rz(wire, angle), kind 1 = cx(c, t)
    let mut ops: Vec<(u8, usize, f64)> = vec![];
    let mut pending: Vec<(u64, f64)> = supports.clone();
    let mut masks = start_masks.clone();
    let mut pi = 0usize;
    while !pending.is_empty() {
        let mut progressed = false;
        let supp: Vec<u64> = pending.iter().map(|(s, _)| *s).collect();
        for s in supp {
            if let Some(w) = (0..n).find(|w| masks[*w] == s) {
                if let Some(pos) = pending.iter().position(|(ms, _)| *ms == s) {
                    let (_, a) = pending[pos];
                    ops.push((0u8, w, a));
                    pending.remove(pos);
                    progressed = true;
                }
            }
        }
        if pending.is_empty() { break; }
        if pi < path.len() {
            let (c, t) = path[pi];
            ops.push((1u8, c, t as f64));
            let mc = masks[c];
            masks[t] ^= mc;
            pi += 1;
            progressed = true;
        }
        if !progressed {
            return Err(PyValueError::new_err("stalled"));
        }
    }
    // uncompute: reverse of the emitted cx prefix
    let cx_prefix: Vec<(usize, usize)> = ops.iter()
        .filter(|(k, _, _)| *k == 1)
        .map(|(_, c, t)| (*c, *t as usize))
        .collect();
    for (c, t) in cx_prefix.iter().rev() {
        ops.push((1u8, *c, *t as f64));
    }
    Ok(ops)
}#[pymodule]
fn compactq_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(block_unitary, m)?)?;
    m.add_function(wrap_pyfunction!(det4, m)?)?;
    m.add_function(wrap_pyfunction!(sim_unitary, m)?)?;
    m.add_function(wrap_pyfunction!(trace2, m)?)?;
    m.add_function(wrap_pyfunction!(py_resynth_1q, m)?)?;
    m.add_function(wrap_pyfunction!(parity_network, m)?)?;
    Ok(())
}
