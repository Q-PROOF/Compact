//! Mirror of compactq.linalg.resynth (and its same_up_to_phase helper).
//!
//! Must track the Python implementation's branch structure: compactq treats the
//! native result as advisory and falls back to Python on any error, and the
//! test suite compares gate decisions against the Python reference.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

#[derive(Clone, Copy)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    fn abs(self) -> f64 {
        self.norm2().sqrt()
    }
    fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }
    fn sub(self, o: Cx) -> Cx {
        Cx {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }
    fn sqrt(self) -> Cx {
        let r = self.norm2().sqrt();
        let th = self.arg() / 2.0;
        Cx {
            re: r * th.cos(),
            im: r * th.sin(),
        }
    }
}

impl std::ops::Add for Cx {
    type Output = Cx;
    fn add(self, o: Cx) -> Cx {
        Cx { re: self.re + o.re, im: self.im + o.im }
    }
}
const fn c(re: f64, im: f64) -> Cx {
    Cx { re, im }
}

fn wrap_angle(x: f64) -> f64 {
    let two_pi = std::f64::consts::TAU;
    let pi = std::f64::consts::PI;
    let mut x = (x + pi) % two_pi;
    if x <= 0.0 {
        x += two_pi;
    }
    x - pi
}

fn same_up_to_phase2(a: [Cx; 4], b: [Cx; 4]) -> bool {
    let mut tr_r = 0.0f64;
    let mut tr_i = 0.0f64;
    for i in 0..4 {
        tr_r += a[i].re * b[i].re + a[i].im * b[i].im;
        tr_i += a[i].im * b[i].re - a[i].re * b[i].im;
    }
    let na: f64 = a.iter().map(|z| z.norm2()).sum();
    let nb: f64 = b.iter().map(|z| z.norm2()).sum();
    let denom = na * nb;
    if denom < 1e-24 {
        return false;
    }
    (tr_r * tr_r + tr_i * tr_i) > (1.0 - 1e-8) * (1.0 - 1e-8) * denom
}

/// Resynthesize a 1-qubit product (4 complex entries, 8 interleaved f64)
/// into the fewest gates; mirrors compactq.linalg.resynth.
#[pyfunction]
#[pyo3(signature = (m, prefer_u))]
pub fn resynth_1q(m: Vec<f64>, prefer_u: bool) -> PyResult<Vec<(String, Vec<f64>)>> {
    if m.len() != 8 {
        return Err(PyValueError::new_err("expected 8 floats"));
    }
    let g = |i: usize| Cx {
        re: m[2 * i],
        im: m[2 * i + 1],
    };
    let prod = [g(0), g(1), g(2), g(3)];
    let sq = 1.0 / 2.0f64.sqrt();
    let named: [(&str, [Cx; 4]); 10] = [
        ("x", [c(0.0, 0.0), c(1.0, 0.0), c(1.0, 0.0), c(0.0, 0.0)]),
        ("y", [c(0.0, 0.0), c(0.0, -1.0), c(0.0, 1.0), c(0.0, 0.0)]),
        ("z", [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(-1.0, 0.0)]),
        ("h", [c(sq, 0.0), c(sq, 0.0), c(sq, 0.0), c(-sq, 0.0)]),
        ("s", [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, 1.0)]),
        ("sdg", [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(0.0, -1.0)]),
        (
            "t",
            [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(sq, sq)],
        ),
        (
            "tdg",
            [c(1.0, 0.0), c(0.0, 0.0), c(0.0, 0.0), c(sq, -sq)],
        ),
        ("sx", [c(0.5, 0.5), c(0.5, -0.5), c(0.5, -0.5), c(0.5, 0.5)]),
        ("sxdg", [c(0.5, -0.5), c(0.5, 0.5), c(0.5, 0.5), c(0.5, -0.5)]),
    ];
    for (name, nm) in named.iter() {
        if same_up_to_phase2(prod, *nm) {
            return Ok(vec![(name.to_string(), vec![])]);
        }
    }

    let (va, vb, vc, vd) = (prod[0], prod[1], prod[2], prod[3]);
    let det = Cx {
        re: va.re * vd.re - va.im * vd.im - (vb.re * vc.re - vb.im * vc.im),
        im: va.re * vd.im + va.im * vd.re - (vb.re * vc.im + vb.im * vc.re),
    };
    if det.abs() < 1e-15 {
        return Err(PyValueError::new_err("non-unitary product"));
    }
    let sroot = det.sqrt();
    let denom = sroot.re * sroot.re + sroot.im * sroot.im;
    let div = |z: Cx| Cx {
        re: (z.re * sroot.re + z.im * sroot.im) / denom,
        im: (z.im * sroot.re - z.re * sroot.im) / denom,
    };
    let va = div(va);
    let vb = div(vb);
    let vc = div(vc);
    let vd = div(vd);

    let tol = 1e-9f64;
    // single-gate RX recognition: V = [[c, -i s], [-i s, c]]
    if (va.re - vd.re).abs() < tol
        && (vb.re - vc.re).abs() < tol
        && vb.re.abs() < tol
        && va.im.abs() < tol
    {
        let theta = wrap_angle(2.0 * (-vb.im).atan2(va.re));
        if theta.abs() < 1e-12 {
            return Ok(vec![]);
        }
        return Ok(vec![("rx".to_string(), vec![theta])]);
    }
    // single-gate RY recognition: V real, [[c, -s], [s, c]]
    if va.im.abs() < tol && vb.im.abs() < tol && vc.im.abs() < tol && vd.im.abs() < tol {
        let theta = wrap_angle(2.0 * vc.re.atan2(va.re));
        if theta.abs() < 1e-12 {
            return Ok(vec![]);
        }
        return Ok(vec![("ry".to_string(), vec![theta])]);
    }

    // exact 2-gate H-times-phase factorisations
    let rsq = 1.0 / 2.0f64.sqrt();
    if (va.abs() - rsq).abs() < 1e-7 && (vb.abs() - rsq).abs() < 1e-7 {
        // H.P(t): V00 = V10, V01 = -V11
        if va.sub(vc).abs() < tol && (vb + vd).abs() < tol {
            let phi = wrap_angle(vb.arg() - va.arg());
            if phi.abs() < 1e-12 {
                return Ok(vec![("h".to_string(), vec![])]);
            }
            return Ok(vec![
                ("p".to_string(), vec![phi]),
                ("h".to_string(), vec![]),
            ]);
        }
        // P.H: V00 = V01, V10 = -V11
        if va.sub(vb).abs() < tol && (vc + vd).abs() < tol {
            let phi = wrap_angle(vc.arg() - va.arg());
            if phi.abs() < 1e-12 {
                return Ok(vec![("h".to_string(), vec![])]);
            }
            return Ok(vec![
                ("h".to_string(), vec![]),
                ("p".to_string(), vec![phi]),
            ]);
        }
    }

    let cos2 = 1.0f64.min(va.abs());
    let theta = 2.0 * vb.abs().atan2(cos2);

    if theta.abs() < 1e-12 {
        let ang = wrap_angle(-2.0 * va.arg());
        if ang.abs() > 1e-12 {
            return Ok(vec![("p".to_string(), vec![ang])]);
        }
        return Ok(vec![]);
    }

    // NB: no 2*pi wrapping on lam_mu/lam_mu_p - RZ is only 4*pi-periodic in
    // each angle and wrapping individually would change the unitary (this is
    // documented in the Python reference and was a real bug in this port).
    let lam_mu = -2.0 * Cx { re: -vb.re, im: -vb.im }.arg();
    let (lam, mu) = if cos2.abs() > 1e-12 {
        let lam_mu_p = -2.0 * va.arg();
        ((lam_mu + lam_mu_p) / 2.0, (lam_mu_p - lam_mu) / 2.0)
    } else {
        (lam_mu, 0.0)
    };

    if prefer_u {
        return Ok(vec![("u3".to_string(), vec![theta, lam, mu])]);
    }
    let mut out: Vec<(String, Vec<f64>)> = Vec::new();
    if mu.abs() > 1e-12 {
        out.push(("rz".to_string(), vec![mu]));
    }
    out.push(("ry".to_string(), vec![theta]));
    if lam.abs() > 1e-12 {
        out.push(("rz".to_string(), vec![lam]));
    }
    Ok(out)
}
