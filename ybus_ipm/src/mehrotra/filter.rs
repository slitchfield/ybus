//! The line search's pieces (`LineSearch::MonotoneFilter`, `docs/LineSearchChecklist.md`):
//! Wächter and Biegler's monotone barrier, from "On the implementation of an interior-point
//! filter line-search algorithm for large-scale nonlinear programming," *Math. Program.*
//! 106(1):25–57, 2006, doi:10.1007/s10107-004-0559-y, with the paper's defaults.
//!
//! - **The barrier problem** at a fixed `μ`: minimize `φ_μ = f − μ (Σ ln z + Σ ln s_l +
//!   Σ ln s_u)` subject to the constraints, `[g; h + z] = 0`.
//! - **Its update** (§3.1): once the barrier problem's error is at most `κ_ε μ`,
//!   `μ ← max(μ_min, min(κ_μ μ, μ^θ_μ))`, as often as the error allows.
//! - **The fraction to the boundary** (§2.2): `τ = max(τ_min, 1 − μ)`.
//! - **The multipliers' safeguard** (§2.2, eq. 16): after each step, every inequality and bound
//!   multiplier is kept within `[μ / (κ_Σ s), κ_Σ μ / s]` of its slack `s`.
//! - **The filter line search** (§2.3): a trial point is accepted if the filter, a set of
//!   `(θ, φ)` pairs, doesn't dominate it, and it either satisfies the Armijo condition on `φ`,
//!   where the switching condition says the step is a descent step for `φ` and `θ` is small, or
//!   reduces `θ` or `φ` enough. Otherwise the step is halved, down to `α_min`. A step that wasn't
//!   an Armijo step adds the current point's margin to the filter.

/// The first barrier parameter, `μ₀`
pub(super) const MU_INIT: f64 = 0.1;
/// How small the barrier problem's error must be, relative to `μ`, before `μ` falls: `κ_ε`
pub(super) const KAPPA_EPSILON: f64 = 10.0;
/// `μ`'s linear and superlinear decreases: `κ_μ` and `θ_μ`
pub(super) const KAPPA_MU: f64 = 0.2;
pub(super) const THETA_MU: f64 = 1.5;
/// The fraction to the boundary's smallest value, `τ_min`
pub(super) const TAU_MIN: f64 = 0.99;
/// How far a multiplier may drift from its slack's `μ / s`: `κ_Σ`
pub(super) const KAPPA_SIGMA: f64 = 1e10;

/// The filter's margins on `θ` and `φ`: `γ_θ` and `γ_φ`
pub(super) const GAMMA_THETA: f64 = 1e-5;
pub(super) const GAMMA_PHI: f64 = 1e-8;
/// The switching condition's constants: `δ`, `s_θ` and `s_φ`
pub(super) const DELTA: f64 = 1.0;
pub(super) const S_THETA: f64 = 1.1;
pub(super) const S_PHI: f64 = 2.3;
/// The Armijo condition's constant, `η_φ`
pub(super) const ETA_PHI: f64 = 1e-4;
/// The safety factor on the smallest step length, `γ_α`
pub(super) const GAMMA_ALPHA: f64 = 0.05;
/// The second-order correction's required decrease in `θ`, `κ_soc`, and most corrections,
///   `p_max`
pub(super) const KAPPA_SOC: f64 = 0.99;
pub(super) const P_MAX: usize = 4;

/// The filter: the `(θ, φ)` pairs a trial point must not be dominated by, with the bound
///   `θ < θ_max` every trial point must meet
pub(super) struct Filter {
    theta_max: f64,
    entries: Vec<(f64, f64)>,
}

impl Filter {
    pub(super) fn new(theta_max: f64) -> Filter {
        Filter {
            theta_max,
            entries: Vec::new(),
        }
    }

    /// Back to its start, `{(θ, φ) : θ ≥ θ_max}`, as when the barrier parameter changes
    pub(super) fn reset(&mut self) {
        self.entries.clear();
    }

    /// Whether `(θ, φ)` is acceptable: below `θ_max`, and with a smaller `θ` or `φ` than each
    ///   entry
    pub(super) fn acceptable(&self, theta: f64, phi: f64) -> bool {
        theta < self.theta_max && self.entries.iter().all(|&(t, f)| theta < t || phi < f)
    }

    /// Adds the margin around the current point, `((1 − γ_θ) θ, φ − γ_φ θ)`, dropping the
    ///   entries it dominates
    pub(super) fn augment(&mut self, theta: f64, phi: f64) {
        let entry = ((1.0 - GAMMA_THETA) * theta, phi - GAMMA_PHI * theta);
        self.entries
            .retain(|&(t, f)| !(t >= entry.0 && f >= entry.1));
        self.entries.push(entry);
    }

    /// Its entries
    pub(super) fn len(&self) -> usize {
        self.entries.len()
    }
}

/// The switching condition (eq. 19): the step is a descent direction for `φ`, by enough
///   relative to `θ`, at step length `alpha`
pub(super) fn switching(alpha: f64, grad_phi_d: f64, theta: f64) -> bool {
    grad_phi_d < 0.0 && alpha * (-grad_phi_d).powf(S_PHI) > DELTA * theta.powf(S_THETA)
}

/// The Armijo condition on `φ` (eq. 20)
pub(super) fn armijo(phi_trial: f64, phi: f64, alpha: f64, grad_phi_d: f64) -> bool {
    phi_trial <= phi + ETA_PHI * alpha * grad_phi_d
}

/// Sufficient decrease of `θ` or `φ` (eq. 18)
pub(super) fn sufficient_decrease(theta_trial: f64, phi_trial: f64, theta: f64, phi: f64) -> bool {
    theta_trial <= (1.0 - GAMMA_THETA) * theta || phi_trial <= phi - GAMMA_PHI * theta
}

/// The smallest step length the line search tries (eq. 23), below which restoration takes over
pub(super) fn alpha_min(theta: f64, grad_phi_d: f64, theta_min: f64) -> f64 {
    let bound = if grad_phi_d < 0.0 && theta <= theta_min {
        GAMMA_THETA
            .min(GAMMA_PHI * theta / -grad_phi_d)
            .min(DELTA * theta.powf(S_THETA) / (-grad_phi_d).powf(S_PHI))
    } else if grad_phi_d < 0.0 {
        GAMMA_THETA.min(GAMMA_PHI * theta / -grad_phi_d)
    } else {
        GAMMA_THETA
    };
    GAMMA_ALPHA * bound
}

/// The barrier function `φ_μ = f − μ (Σ ln z + Σ ln s_l + Σ ln s_u)`, for the objective `f`
///   and the slacks given
pub(super) fn barrier(f: f64, mu: f64, slacks: [&[f64]; 3]) -> f64 {
    f - mu
        * slacks
            .iter()
            .flat_map(|s| s.iter())
            .map(|s| s.ln())
            .sum::<f64>()
}

/// The next barrier parameter after the barrier problem's error fell to at most `κ_ε μ`
pub(super) fn decrease(mu: f64, floor: f64) -> f64 {
    (KAPPA_MU * mu).min(mu.powf(THETA_MU)).max(floor)
}

/// The fraction to the boundary at barrier parameter `mu`
pub(super) fn fraction(mu: f64) -> f64 {
    (1.0 - mu).max(TAU_MIN)
}

/// Keeps each multiplier in `multipliers` within `[μ / (κ_Σ s), κ_Σ μ / s]` of its slack in
///   `slacks` (§2.2, eq. 16), so that the primal-dual Hessian stays near the primal one
pub(super) fn safeguard(multipliers: &mut [f64], slacks: &[f64], mu: f64) {
    for (v, &s) in multipliers.iter_mut().zip(slacks) {
        *v = v.clamp(mu / (KAPPA_SIGMA * s), KAPPA_SIGMA * mu / s);
    }
}

/// The barrier problem's complementarity error, `max |s v − μ|` over the slacks and
///   multipliers given
pub(super) fn complementarity_error(pairs: [(&[f64], &[f64]); 3], mu: f64) -> f64 {
    pairs
        .iter()
        .flat_map(|(s, v)| s.iter().zip(v.iter()))
        .map(|(s, v)| (s * v - mu).abs())
        .fold(0.0, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_barrier_parameter_falls_linearly_then_superlinearly_to_its_floor() {
        // κ_μ μ while it's the smaller, μ^1.5 once that is, and never below the floor
        assert!((decrease(0.1, 0.0) - 0.02).abs() < 1e-15);
        assert!((decrease(1e-4, 0.0) - 1e-6).abs() < 1e-18);
        assert_eq!(decrease(1e-4, 1e-5), 1e-5);
    }

    #[test]
    fn the_filter_rejects_what_its_entries_dominate_and_keeps_only_undominated_ones() {
        let mut f = Filter::new(10.0);
        assert!(f.acceptable(5.0, 100.0));
        assert!(!f.acceptable(10.0, -1e10), "θ must be below θ_max");
        f.augment(4.0, 3.0);
        let (t, p) = ((1.0 - GAMMA_THETA) * 4.0, 3.0 - GAMMA_PHI * 4.0);
        assert!(!f.acceptable(t, p), "the entry itself");
        assert!(!f.acceptable(5.0, 4.0), "worse in both");
        assert!(f.acceptable(1.0, 100.0), "better θ");
        assert!(f.acceptable(9.0, -1.0), "better φ");
        // A dominating margin replaces the entries it dominates
        f.augment(5.0, 10.0);
        assert_eq!(f.len(), 2);
        f.augment(1.0, 1.0);
        assert_eq!(f.len(), 1, "the margin around (1, 1) dominates both");
        f.reset();
        assert!(f.acceptable(5.0, 4.0));
    }

    #[test]
    fn the_acceptance_tests() {
        // Descent, by enough against θ, switches; not at a small step or an ascent
        assert!(switching(1.0, -1.0, 0.5));
        assert!(!switching(1e-3, -1.0, 0.5));
        assert!(!switching(1.0, 1.0, 0.0));
        assert!(armijo(1.0 - 2e-4, 1.0, 1.0, -1.0));
        assert!(!armijo(1.0, 1.0, 1.0, -1.0));
        assert!(sufficient_decrease(0.5, 9.0, 1.0, 1.0), "θ halved");
        assert!(sufficient_decrease(1.0, 0.5, 1.0, 1.0), "φ down");
        assert!(!sufficient_decrease(1.0, 1.0, 1.0, 1.0));
        // α_min is γ_α γ_θ without descent, and smaller with it
        assert_eq!(alpha_min(1.0, 1.0, 1e-4), GAMMA_ALPHA * GAMMA_THETA);
        assert!(alpha_min(1e-5, -1.0, 1e-4) < GAMMA_ALPHA * GAMMA_THETA);
    }

    #[test]
    fn the_fraction_to_the_boundary_nears_one_as_the_barrier_falls() {
        assert_eq!(fraction(0.1), 0.99);
        assert_eq!(fraction(1e-4), 1.0 - 1e-4);
    }

    #[test]
    fn the_safeguard_clamps_only_multipliers_far_from_their_slacks() {
        let mut v = [1.0, 1e-20, 1e20];
        safeguard(&mut v, &[1.0, 1.0, 1.0], 0.1);
        assert_eq!(v, [1.0, 0.1 / KAPPA_SIGMA, KAPPA_SIGMA * 0.1]);
    }

    #[test]
    fn the_barrier_function_and_complementarity_error() {
        let one = [1.0];
        let e = std::f64::consts::E;
        assert!((barrier(2.0, 0.5, [&[e], &one, &[]]) - 1.5).abs() < 1e-15);
        assert_eq!(
            complementarity_error([(&[2.0], &[3.0]), (&[], &[]), (&[], &[])], 1.0),
            5.0
        );
    }
}
