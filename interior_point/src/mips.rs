//! MIPS, the MATPOWER Interior Point Solver: a primal-dual interior-point method for nonlinear
//! programs, reimplemented from MATPOWER's BSD-licensed `mips.m` (MIPS 1.5.x, in MATPOWER 8.1)
//! so that its iterations, results and multipliers match MATPOWER's (Wang, Murillo-Sánchez,
//! Zimmerman and Thomas, 2007).

/// MIPS's options, with its defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct MipsOptions {
    /// Termination tolerance for feasibility (`feastol`)
    pub feastol: f64,
    /// Termination tolerance for the gradient (`gradtol`)
    pub gradtol: f64,
    /// Termination tolerance for complementarity (`comptol`)
    pub comptol: f64,
    /// Termination tolerance for the cost (`costtol`)
    pub costtol: f64,
    /// The most iterations (`max_it`)
    pub max_it: usize,
    /// Whether to reduce a step that worsens both feasibility and the gradient
    /// (`step_control`)
    pub step_control: bool,
    /// The most step-size reductions per step, with step control (`sc.red_it`)
    pub red_it: usize,
    /// The objective's scale factor (`cost_mult`)
    pub cost_mult: f64,
    /// The fraction of the distance to the boundary a step may take (`xi`), in `[0.5, 1)`
    pub xi: f64,
    /// The centering parameter (`sigma`), in `(0, 1]`
    pub sigma: f64,
    /// The slack variables' starting value (`z0`)
    pub z0: f64,
    /// The smallest acceptable primal or dual step length (`alpha_min`)
    pub alpha_min: f64,
    /// The bounds on step control's ratio of actual to predicted change (`rho_min`, `rho_max`)
    pub rho_min: f64,
    /// See `rho_min`
    pub rho_max: f64,
    /// The threshold below which an inequality's multiplier is set to zero when the inequality
    /// isn't binding (`mu_threshold`)
    pub mu_threshold: f64,
    /// The largest acceptable 2-norm of the reduced Newton step (`max_stepsize`)
    pub max_stepsize: f64,
}

impl Default for MipsOptions {
    fn default() -> MipsOptions {
        MipsOptions {
            feastol: 1e-6,
            gradtol: 1e-6,
            comptol: 1e-6,
            costtol: 1e-6,
            max_it: 150,
            step_control: false,
            red_it: 20,
            cost_mult: 1.0,
            xi: 0.99995,
            sigma: 0.1,
            z0: 1.0,
            alpha_min: 1e-8,
            rho_min: 0.95,
            rho_max: 1.05,
            mu_threshold: 1e-5,
            max_stepsize: 1e10,
        }
    }
}

/// How MIPS ended (`exitflag`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MipsStatus {
    /// Every termination test passed (`exitflag` 1)
    Converged,
    /// `max_it` iterations passed without converging (0)
    MaxIterations,
    /// A step was NaN or too long, a step length or the barrier parameter left its range
    /// (−1)
    NumericalFailure,
}

/// The multipliers at the solution, as MIPS reports them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Multipliers {
    /// One per nonlinear equality constraint (`eqnonlin`)
    pub eq_nonlin: Vec<f64>,
    /// One per nonlinear inequality constraint (`ineqnonlin`)
    pub ineq_nonlin: Vec<f64>,
    /// One per linear constraint, for its lower bound `l` (`mu_l`)
    pub mu_l: Vec<f64>,
    /// One per linear constraint, for its upper bound `u` (`mu_u`)
    pub mu_u: Vec<f64>,
    /// One per variable, for `xmin` (`lower`)
    pub lower: Vec<f64>,
    /// One per variable, for `xmax` (`upper`)
    pub upper: Vec<f64>,
}

/// One iteration's record (`output.hist`): the termination measures, the barrier parameter,
/// the step's 2-norm, the objective and the step lengths.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IterationRecord {
    /// The feasibility measure
    pub feascond: f64,
    /// The gradient measure
    pub gradcond: f64,
    /// The complementarity measure
    pub compcond: f64,
    /// The cost measure
    pub costcond: f64,
    /// The barrier parameter, `γ`
    pub gamma: f64,
    /// The Newton step's 2-norm, `‖dx‖`
    pub stepsize: f64,
    /// The objective, without `cost_mult`
    pub obj: f64,
    /// The primal step length
    pub alphap: f64,
    /// The dual step length
    pub alphad: f64,
}

/// MIPS's result.
#[derive(Clone, Debug, PartialEq)]
pub struct MipsSolution {
    /// The final point
    pub x: Vec<f64>,
    /// The objective there
    pub f: f64,
    /// How the solve ended
    pub status: MipsStatus,
    /// The number of iterations taken
    pub iterations: usize,
    /// The multipliers
    pub multipliers: Multipliers,
    /// One record for the starting point and one per iteration
    pub history: Vec<IterationRecord>,
}

impl MipsSolution {
    /// Whether the solve converged.
    pub fn converged(&self) -> bool {
        self.status == MipsStatus::Converged
    }
}
