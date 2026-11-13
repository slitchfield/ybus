# ybus_ipm

Interior-point solvers for nonlinear programs:

minimize `f(x)` subject to `g(x) = 0`, `h(x) ≤ 0`, `l ≤ A x ≤ u` and `xmin ≤ x ≤ xmax`.

It has two solvers, which take the same problems and return the same `Solution`:

- **MIPS**, the MATPOWER Interior Point Solver, reimplemented step for step from MATPOWER's
  `mips.m` (`mips::mips`, and `mips::qps_mips` for quadratic and linear programs). Its
  iterations and results are checked against MATPOWER's own, problem by problem.
- **A primal-dual method with Mehrotra's predictor-corrector**, written from the literature
  (`mehrotra::mehrotra` and `mehrotra::qps_mehrotra`). It usually takes fewer iterations.

A general problem implements the `Problem` trait, which supplies the objective and any
nonlinear constraints with their first and second derivatives. The KKT systems are factored by
LU or by LDLᵀ, using the sparse matrices of
[`ybus_sparse`](https://crates.io/crates/ybus_sparse).

It is part of the [`ybus`](https://github.com/slitchfield/ybus) workspace, where it solves the
DC and AC optimal power flows of
[`ybus_powerflow`](https://crates.io/crates/ybus_powerflow).

## Status

0.1, the first release. The solvers are tested against MATPOWER's results, but the API is young:
expect breaking changes in the 0.x releases, each one a minor version bump.

## Quick start

```toml
[dependencies]
ybus_ipm = "0.1"
ybus_sparse = "0.1"
```

A quadratic program: minimize `½ xᵀ H x + cᵀ x` subject to three linear constraints and
`x ≥ 0`.

```rust
use ybus_ipm::mips::{qps_mips, MipsOptions};
use ybus_ipm::{Constraints, InteriorPointError, QuadraticProblem};
use ybus_sparse::{CsrMatrix, SparseMatrix};

fn dense(rows: &[&[f64]]) -> CsrMatrix {
    let triplets: Vec<_> = rows
        .iter()
        .enumerate()
        .flat_map(|(i, row)| row.iter().enumerate().map(move |(j, &v)| (i, j, v)))
        .collect();
    SparseMatrix::from_triplets(rows.len(), rows[0].len(), triplets)
        .unwrap()
        .to_csr()
}

fn main() -> Result<(), InteriorPointError> {
    let qp = QuadraticProblem {
        h: Some(dense(&[&[1.0, -1.0], &[-1.0, 2.0]])),
        c: vec![-2.0, -6.0],
    };
    let constraints = Constraints {
        a: Some(dense(&[&[1.0, 1.0], &[-1.0, 2.0], &[2.0, 1.0]])),
        l: None, // A x has no lower bound
        u: Some(vec![2.0, 2.0, 3.0]),
        xmin: Some(vec![0.0; 2]),
        xmax: None,
    };
    let solution = qps_mips(&qp, &constraints, None, &MipsOptions::default())?;
    assert!(solution.converged());
    assert!((solution.x[0] - 2.0 / 3.0).abs() < 1e-6);
    assert!((solution.x[1] - 4.0 / 3.0).abs() < 1e-6);
    // The first two constraints are active, so their multipliers are positive
    assert!(solution.multipliers.mu_u[0] > 0.0 && solution.multipliers.mu_u[1] > 0.0);
    Ok(())
}
```

`mehrotra::qps_mehrotra(&qp, &constraints, None, &MehrotraOptions::default())` solves the same
program with the predictor-corrector method.

## Development

- **MSRV:** Rust 1.88, set for the whole workspace in the root `Cargo.toml`.
- **Tests:** in the repository's `ybus_ipm/tests/`. They compare against MIPS's own test
  problems and MATPOWER's answers, and check the KKT conditions of each result.
- **API documentation:** `cargo doc -p ybus_ipm --open`.

See the [workspace README](https://github.com/slitchfield/ybus#readme) for `./check.sh`, CI and
the rest of the docs, and
[`docs/OptimizationChecklist.md`](https://github.com/slitchfield/ybus/blob/main/docs/OptimizationChecklist.md)
for the solvers' roadmap.

## License

The crate is licensed under either of the Apache License, Version 2.0
([`LICENSE-APACHE`](https://github.com/slitchfield/ybus/blob/main/LICENSE-APACHE)) or the MIT
license ([`LICENSE-MIT`](https://github.com/slitchfield/ybus/blob/main/LICENSE-MIT)), at your
option.

The MIPS port in `src/mips.rs` is also derived from MIPS, and redistributed under MIPS's
3-clause BSD license; see
[`LICENSE-MIPS`](https://github.com/slitchfield/ybus/blob/main/LICENSE-MIPS) for the files it
covers. The crate's license is therefore `(MIT OR Apache-2.0) AND BSD-3-Clause`.

See [`AI-DISCLAIMER.md`](https://github.com/slitchfield/ybus/blob/main/AI-DISCLAIMER.md) for
how this licensing relates to the use of AI in developing the crate.
