//! The reduced KKT system both solvers solve at each iteration,
//!
//! `[W Jgᵀ; Jg 0] [dx; dλ] = [r_x; r_λ]`,
//!
//! factored by LU, the default, or by a pivoted LDLᵀ ([`KktFactorization`]).
//!
//! **LU** ([`KKT_ORDERING`]) factors the matrix as assembled, as
//! MATPOWER's MIPS does with its backslash. MIPS's results match MATPOWER's with it.
//!
//! **LDLᵀ** is `ybus_sparse`'s multifrontal factorization with threshold Bunch–Kaufman pivoting
//! ([`Ldlt`]), ordered by AMD:
//!
//! - **Exactly symmetric:** the lower triangle of `W` is mirrored, so the factorization's symmetry
//!   check holds, and every diagonal is stored, so the pattern stays fixed and the analysis is
//!   kept.
//! - **No regularization:** 2×2 pivots handle the zero diagonals, of the constraints' block and of
//!   `W`'s, that defeat static pivots.
//! - **The threshold, 1e-4:** measured on the KKT matrices of DC and AC `case9241pegase` and
//!   `case_ACTIVSg25k`. MA57's 0.01 delayed so many pivots that the factor outgrew LU's, and 1e-8
//!   lost accuracy that refinement couldn't recover.
//! - **Refined:** each solve takes up to 3 steps of normwise refinement against the matrix, until
//!   the relative residual is at most 1e-12 or stops halving; it reached LU's residual on those
//!   matrices.
//! - **The inertia:** the factorization gives it, and the predictor-corrector method corrects an
//!   indefinite `W` with it (see [`crate::mehrotra`]). MIPS only swaps the factorization.

use std::cell::Cell;
use std::rc::Rc;

use ybus_sparse::compressed_matrix::{
    AssemblyMap, Csc, CscMatrix, CsrMatrix, Inertia, Ldlt, LdltOptions, Lu, LuOptions, Ordering,
    SymbolicLdlt, SymbolicLu,
};
use ybus_sparse::sparse_error::SparseError;
use ybus_sparse::Factorization;

use crate::mips::KKT_ORDERING;
use crate::newton::Pattern;

/// How the solvers factor the KKT matrix (see the [module documentation](self)).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum KktFactorization {
    /// LU, ordered by [`KKT_ORDERING`], as assembled: MATPOWER's
    /// choice, and the default
    #[default]
    Lu,
    /// The pivoted LDLᵀ, ordered by approximate minimum degree, with each solve refined. It
    /// reports the matrix's inertia, which the predictor-corrector method uses to correct an
    /// indefinite Hessian.
    Ldlt,
}

/// LDLᵀ's pivot threshold
const THRESHOLD: f64 = 1e-4;
/// LDLᵀ's refinement: the most steps, and the relative residual it stops at
const REFINEMENT_STEPS: usize = 3;
const REFINEMENT_TOLERANCE: f64 = 1e-12;

/// The KKT matrix's factorization, with its analysis kept across a solve's iterations while
/// the pattern holds, and counts of the analyses and refinement steps.
pub(crate) struct KktSolver {
    method: KktFactorization,
    lu: Option<SymbolicLu>,
    ldlt: Option<SymbolicLdlt>,
    analyses: usize,
    /// Shared with each factor, which adds its solves' steps
    refinement_steps: Rc<Cell<usize>>,
    /// The KKT matrix's assembly, kept while its inputs' patterns hold
    assembly: Option<Assembly>,
}

/// The KKT matrix's assembly for one set of input patterns: `from_triplets`' sums replayed
/// by kept [`AssemblyMap`]s, bit for bit (see `KktSolver::factor`)
struct Assembly {
    nx: usize,
    w: Pattern,
    bounds: Vec<usize>,
    jg: Pattern,
    kind: AssemblyKind,
}

enum AssemblyKind {
    /// LU's CSC matrix from `W`'s and the bounds' triplets, then `Jg`'s twice
    Lu {
        matrix: CscMatrix<f64>,
        map: AssemblyMap<Csc>,
    },
    /// LDLᵀ's lower triangle, then the matrix with it mirrored
    Ldlt(Box<LdltAssembly>),
}

struct LdltAssembly {
    /// The slots of `W`'s values, then the bounds' (numbered after them), kept in the lower
    /// triangle, in order
    lower_sources: Vec<usize>,
    lower: CsrMatrix<f64>,
    lower_map: AssemblyMap,
    /// The full matrix's triplets' sources in the lower triangle's values
    mirror_sources: Vec<usize>,
    matrix: CsrMatrix<f64>,
    map: AssemblyMap,
}

impl Assembly {
    fn matches(
        &self,
        nx: usize,
        w: &CsrMatrix<f64>,
        bounds: &[(usize, f64)],
        jg: &CsrMatrix<f64>,
    ) -> bool {
        self.nx == nx
            && self.w.matches(w)
            && self.jg.matches(jg)
            && self.bounds.len() == bounds.len()
            && self.bounds.iter().zip(bounds).all(|(&k, &(b, _))| k == b)
    }
}

/// An assembled KKT matrix, in the orientation its factorization takes
enum Assembled {
    Lu(CscMatrix<f64>),
    Ldlt(CsrMatrix<f64>),
}

/// `W`'s triplets in row order, then the bounds' terms on the diagonal: the order the solvers
/// have always assembled them in
fn w_triplets<'a>(
    w: &'a CsrMatrix<f64>,
    bounds: &'a [(usize, f64)],
) -> impl Iterator<Item = (usize, usize, f64)> + 'a {
    w.iter().chain(bounds.iter().map(|&(k, v)| (k, k, v)))
}

/// One iteration's factored KKT matrix.
pub(crate) struct KktFactor {
    kind: FactorKind,
    refinement_steps: Rc<Cell<usize>>,
}

enum FactorKind {
    Lu(Box<Lu<f64>>),
    Ldlt {
        ldlt: Ldlt<f64>,
        /// The factored matrix, which the solves refine against
        matrix: CsrMatrix<f64>,
    },
}

impl KktSolver {
    pub(crate) fn new(method: KktFactorization) -> KktSolver {
        KktSolver {
            method,
            lu: None,
            ldlt: None,
            analyses: 0,
            refinement_steps: Rc::new(Cell::new(0)),
            assembly: None,
        }
    }

    /// Factors `[W Jgᵀ; Jg 0]` (`nx + neq` square), from `W` (`nx × nx`) plus the bounds' terms
    /// `bounds` on its diagonal, `(k, value)` added at `(k, k)` after `W`'s entries, and `Jg`
    /// (`neq × nx`), with `delta_w` added to `W`'s diagonal for LDLᵀ (an inertia correction; LU
    /// ignores it). Returns the analysis's and factorization's errors.
    ///
    /// The matrix is assembled as the solvers always have, from triplets with duplicates
    /// summed: `W`'s in row order and the bounds', then each entry of `Jg` at `(j, nx + k)` and
    /// `(nx + k, j)` for LU; for LDLᵀ, those of them in the lower triangle, `Jg`'s at
    /// `(nx + k, j)`, and every diagonal, then mirrored. The sums are replayed by kept
    /// assembly maps while the inputs' patterns hold, bit for bit what `from_triplets` gives.
    pub(crate) fn factor(
        &mut self,
        nx: usize,
        w: &CsrMatrix<f64>,
        bounds: &[(usize, f64)],
        jg: &CsrMatrix<f64>,
        delta_w: f64,
    ) -> Result<KktFactor, SparseError> {
        let kind = match self.assemble(nx, w, bounds, jg, delta_w) {
            Assembled::Lu(kkt) => FactorKind::Lu(Box::new(self.factor_lu(&kkt)?)),
            Assembled::Ldlt(matrix) => self.factor_ldlt(matrix)?,
        };
        Ok(KktFactor {
            kind,
            refinement_steps: Rc::clone(&self.refinement_steps),
        })
    }

    /// The KKT matrix `factor` factors (see there), from the kept assembly maps, planned afresh
    /// if the inputs' patterns differ from the last ones
    fn assemble(
        &mut self,
        nx: usize,
        w: &CsrMatrix<f64>,
        bounds: &[(usize, f64)],
        jg: &CsrMatrix<f64>,
        delta_w: f64,
    ) -> Assembled {
        if !self
            .assembly
            .as_ref()
            .is_some_and(|a| a.matches(nx, w, bounds, jg))
        {
            self.assembly = Some(self.plan_assembly(nx, w, bounds, jg));
        }
        match &self.assembly.as_ref().expect("just planned").kind {
            AssemblyKind::Lu { matrix, map } => {
                let mut values: Vec<f64> = w_triplets(w, bounds).map(|(_, _, v)| v).collect();
                for &v in jg.values() {
                    values.extend([v, v]);
                }
                let mut kkt = matrix.clone();
                kkt.refill(map, &values).expect("planned");
                Assembled::Lu(kkt)
            }
            AssemblyKind::Ldlt(assembly) => {
                let LdltAssembly {
                    lower_sources,
                    lower,
                    lower_map,
                    mirror_sources,
                    matrix,
                    map,
                } = assembly.as_ref();
                let size = nx + jg.shape().0;
                let wb: Vec<f64> = w_triplets(w, bounds).map(|(_, _, v)| v).collect();
                let mut values: Vec<f64> = lower_sources.iter().map(|&s| wb[s]).collect();
                values.extend_from_slice(jg.values());
                values.extend((0..size).map(|i| if i < nx { delta_w } else { 0.0 }));
                let mut lower = lower.clone();
                lower.refill(lower_map, &values).expect("planned");
                let both: Vec<f64> = mirror_sources.iter().map(|&s| lower.values()[s]).collect();
                let mut matrix = matrix.clone();
                matrix.refill(map, &both).expect("planned");
                Assembled::Ldlt(matrix)
            }
        }
    }

    /// The assembly maps for these inputs' patterns, from their triplets' positions
    fn plan_assembly(
        &self,
        nx: usize,
        w: &CsrMatrix<f64>,
        bounds: &[(usize, f64)],
        jg: &CsrMatrix<f64>,
    ) -> Assembly {
        let size = nx + jg.shape().0;
        let kind = match self.method {
            KktFactorization::Lu => {
                let mut triplets: Vec<(usize, usize, f64)> = w_triplets(w, bounds).collect();
                for (k, j, v) in jg.iter() {
                    triplets.push((j, nx + k, v));
                    triplets.push((nx + k, j, v));
                }
                let (matrix, map) = CscMatrix::from_triplets_mapped(size, size, triplets)
                    .expect("indices in range");
                AssemblyKind::Lu { matrix, map }
            }
            KktFactorization::Ldlt => {
                // The lower triangle, every diagonal stored (with δ_w on W's), then mirrored:
                //   exactly symmetric
                let (lower_sources, mut lower): (Vec<usize>, Vec<(usize, usize, f64)>) =
                    w_triplets(w, bounds)
                        .enumerate()
                        .filter(|&(_, (r, c, _))| r >= c)
                        .unzip();
                lower.extend(jg.iter().map(|(k, j, v)| (nx + k, j, v)));
                lower.extend((0..size).map(|i| (i, i, 0.0)));
                let (lower, lower_map) =
                    CsrMatrix::from_triplets_mapped(size, size, lower).expect("indices in range");
                let mut mirror_sources = Vec::with_capacity(2 * lower.num_stored());
                let mut both = Vec::with_capacity(2 * lower.num_stored());
                for (s, (r, c, v)) in lower.iter().enumerate() {
                    mirror_sources.push(s);
                    both.push((r, c, v));
                    if r != c {
                        mirror_sources.push(s);
                        both.push((c, r, v));
                    }
                }
                let (matrix, map) =
                    CsrMatrix::from_triplets_mapped(size, size, both).expect("indices in range");
                AssemblyKind::Ldlt(Box::new(LdltAssembly {
                    lower_sources,
                    lower,
                    lower_map,
                    mirror_sources,
                    matrix,
                    map,
                }))
            }
        };
        Assembly {
            nx,
            w: Pattern::of(w),
            bounds: bounds.iter().map(|&(k, _)| k).collect(),
            jg: Pattern::of(jg),
            kind,
        }
    }

    /// Factors `kkt` with the kept analysis, or, if there is none yet or `kkt`'s pattern
    /// differs from the one it was made for, with a fresh one, which is then kept.
    fn factor_lu(&mut self, kkt: &CscMatrix<f64>) -> Result<Lu<f64>, SparseError> {
        if let Some(symbolic) = &self.lu {
            match symbolic.factor(kkt, LuOptions::default()) {
                Err(SparseError::PatternMismatch) => {}
                factored => return factored,
            }
        }
        let symbolic = kkt.analyze_lu(KKT_ORDERING)?;
        self.analyses += 1;
        let factored = symbolic.factor(kkt, LuOptions::default());
        self.lu = Some(symbolic);
        factored
    }

    /// Factors the assembled, exactly symmetric `matrix` by LDLᵀ with the kept analysis, or a
    /// fresh one if there is none yet or the pattern differs.
    fn factor_ldlt(&mut self, matrix: CsrMatrix<f64>) -> Result<FactorKind, SparseError> {
        let options = LdltOptions {
            threshold: THRESHOLD,
        };
        if let Some(symbolic) = &self.ldlt {
            match symbolic.factor(&matrix, options) {
                Err(SparseError::PatternMismatch) => {}
                factored => {
                    return Ok(FactorKind::Ldlt {
                        ldlt: factored?,
                        matrix,
                    })
                }
            }
        }
        let symbolic = matrix.analyze_ldlt(Ordering::Amd)?;
        self.analyses += 1;
        let factored = symbolic.factor(&matrix, options);
        self.ldlt = Some(symbolic);
        Ok(FactorKind::Ldlt {
            ldlt: factored?,
            matrix,
        })
    }

    /// The analyses made so far.
    pub(crate) fn analyses(&self) -> usize {
        self.analyses
    }

    /// The refinement steps taken so far, over every solve (0 with LU).
    pub(crate) fn refinement_steps(&self) -> usize {
        self.refinement_steps.get()
    }
}

impl KktFactor {
    /// Solves the KKT system for `rhs`: with LU, directly; with LDLᵀ, refined.
    pub(crate) fn solve(&self, rhs: &[f64]) -> Result<Vec<f64>, SparseError> {
        match &self.kind {
            FactorKind::Lu(lu) => lu.solve(rhs),
            FactorKind::Ldlt { ldlt, matrix } => {
                let refined = ldlt.solve_refined_normwise(
                    matrix,
                    rhs,
                    REFINEMENT_STEPS,
                    REFINEMENT_TOLERANCE,
                )?;
                self.refinement_steps
                    .set(self.refinement_steps.get() + refined.steps);
                Ok(refined.x)
            }
        }
    }

    /// The factored matrix's inertia, with LDLᵀ (`None` with LU).
    pub(crate) fn inertia(&self) -> Option<Inertia> {
        match &self.kind {
            FactorKind::Lu(_) => None,
            FactorKind::Ldlt { ldlt, .. } => Some(ldlt.inertia()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The KKT matrix as the solvers assembled it before the kept maps, straight from triplets
    fn reference(
        method: KktFactorization,
        nx: usize,
        w: &CsrMatrix<f64>,
        bounds: &[(usize, f64)],
        jg: &CsrMatrix<f64>,
        delta_w: f64,
    ) -> (Vec<usize>, Vec<usize>, Vec<u64>) {
        let size = nx + jg.shape().0;
        let w: Vec<(usize, usize, f64)> = w_triplets(w, bounds).collect();
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        match method {
            KktFactorization::Lu => {
                let mut triplets = w;
                for (k, j, v) in jg.iter() {
                    triplets.push((j, nx + k, v));
                    triplets.push((nx + k, j, v));
                }
                let kkt = CscMatrix::from_triplets(size, size, triplets).unwrap();
                (
                    kkt.pointers().to_vec(),
                    kkt.indices().to_vec(),
                    bits(kkt.values()),
                )
            }
            KktFactorization::Ldlt => {
                let mut lower: Vec<(usize, usize, f64)> =
                    w.into_iter().filter(|&(r, c, _)| r >= c).collect();
                lower.extend(jg.iter().map(|(k, j, v)| (nx + k, j, v)));
                lower.extend((0..size).map(|i| (i, i, if i < nx { delta_w } else { 0.0 })));
                let lower = CsrMatrix::from_triplets(size, size, lower).unwrap();
                let mut both = vec![];
                for (r, c, v) in lower.iter() {
                    both.push((r, c, v));
                    if r != c {
                        both.push((c, r, v));
                    }
                }
                let m = CsrMatrix::from_triplets(size, size, both).unwrap();
                (
                    m.pointers().to_vec(),
                    m.indices().to_vec(),
                    bits(m.values()),
                )
            }
        }
    }

    fn assembled(a: Assembled) -> (Vec<usize>, Vec<usize>, Vec<u64>) {
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        match a {
            Assembled::Lu(m) => (
                m.pointers().to_vec(),
                m.indices().to_vec(),
                bits(m.values()),
            ),
            Assembled::Ldlt(m) => (
                m.pointers().to_vec(),
                m.indices().to_vec(),
                bits(m.values()),
            ),
        }
    }

    fn matrix(rows: usize, cols: usize, entries: &[(usize, usize, f64)]) -> CsrMatrix<f64> {
        CsrMatrix::from_triplets(rows, cols, entries.to_vec()).unwrap()
    }

    #[test]
    fn the_kept_assembly_matches_the_triplets_bit_for_bit() {
        // W with a duplicate-prone diagonal (the bounds add to it), Jg with two rows
        let w1 = matrix(
            3,
            3,
            &[
                (0, 0, 2.0),
                (0, 2, -0.5),
                (1, 1, 1e-17),
                (2, 0, -0.5),
                (2, 2, 4.0),
            ],
        );
        let w2 = matrix(
            3,
            3,
            &[
                (0, 0, -0.0),
                (0, 2, 3.0),
                (1, 1, 1e16),
                (2, 0, 3.0),
                (2, 2, 1.0),
            ],
        );
        let jg = matrix(2, 3, &[(0, 0, 1.0), (0, 1, -1.0), (1, 2, 0.25)]);
        // A changed pattern partway, as a problem whose pattern changes gives
        let w3 = matrix(3, 3, &[(0, 0, 1.0), (1, 1, 2.0), (1, 2, 0.5), (2, 1, 0.5)]);
        for method in [KktFactorization::Lu, KktFactorization::Ldlt] {
            let mut solver = KktSolver::new(method);
            for (w, bounds, delta_w) in [
                (&w1, vec![(0, 0.5), (2, 1e-9)], 0.0),
                (&w2, vec![(0, -0.0), (2, 7.0)], 1e-4),
                (&w1, vec![(0, 3.0), (2, -1.0)], 0.0),
                (&w3, vec![(0, 0.5), (2, 1e-9)], 0.0),
                (&w3, vec![(1, 0.5)], 2.0),
                (&w1, vec![], 0.0),
            ] {
                let planned = assembled(solver.assemble(3, w, &bounds, &jg, delta_w));
                assert_eq!(
                    planned,
                    reference(method, 3, w, &bounds, &jg, delta_w),
                    "{method:?}"
                );
            }
        }
    }
}
