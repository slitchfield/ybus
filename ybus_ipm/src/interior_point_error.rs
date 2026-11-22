//! The error type for the solvers' inputs.

use std::fmt;

/// Why a solver couldn't start. A solve that starts and then fails (numerically, or at its
/// iteration limit) is not an error: its result reports how it ended.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum InteriorPointError {
    /// An option is out of its range, such as MIPS's `xi` outside `[0.5, 1)` or `sigma` outside
    /// `(0, 1]`, as `mips.m` checks them.
    InvalidOption {
        /// The option's name
        name: &'static str,
    },
    /// An input's length or shape doesn't match the number of variables or the linear
    /// constraints.
    Dimension {
        /// Which input
        what: &'static str,
    },
    /// A linear program (no objective Hessian) has no constraints or bounds, which MATPOWER's
    /// `qps_mips` refuses.
    UnboundedLinearProgram,
    /// The predictor-corrector method's line search was asked for with the KKT systems factored
    /// by LU: a filter needs LDLᵀ, whose inertia shows whether a step is a descent direction
    /// (`MehrotraOptions::line_search`, `kkt_factorization`).
    LineSearchNeedsLdlt,
}

impl fmt::Display for InteriorPointError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InteriorPointError::InvalidOption { name } => {
                write!(f, "option {name} is out of range")
            }
            InteriorPointError::Dimension { what } => write!(f, "{what} has the wrong dimension"),
            InteriorPointError::UnboundedLinearProgram => {
                write!(
                    f,
                    "a linear program must include constraints or variable bounds"
                )
            }
            InteriorPointError::LineSearchNeedsLdlt => {
                write!(f, "a line search needs the KKT systems factored by LDLᵀ")
            }
        }
    }
}

impl std::error::Error for InteriorPointError {}
