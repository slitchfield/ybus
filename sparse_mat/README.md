# sparse_mat

A sparse matrix library for Rust: construction, arithmetic, products, factorizations and
iterative solvers, Matrix Market I/O, and optional `serde` and `rayon` support. It is generic
over the element type (`f64` by default) and has two required dependencies, `num-traits` and `num-complex`.

It has two matrix types:

- **`SparseMatrix`** builds and edits a matrix: entries live in a hash map, so inserting,
  removing and looking up single entries is cheap, and a zero is never stored.
- **`CompressedMatrix`** does the numerical work: products, solvers, norms, selections and
  arithmetic. It stores its entries compressed along one axis, with the orientation in the
  type: rows in a `CsrMatrix`, columns in a `CscMatrix`. Its pattern (which positions are
  stored) is fixed when it's built, so values can be refilled without changing the structure.

Convert with `to_csr()`/`to_csc()` once the matrix is built, or build a compressed matrix
directly from triplets or raw arrays.

It is the numerical foundation of the [`ybus`](https://github.com/slitchfield/ybus) workspace,
which is building DC and AC power flow on top of it; see
[`docs/PowerFlowChecklist.md`](../docs/PowerFlowChecklist.md) for that roadmap.

## Quick start

```toml
[dependencies]
sparse_mat = { git = "https://github.com/slitchfield/ybus" }
```

```rust
use sparse_mat::compressed_matrix::CsrMatrix;
use sparse_mat::sparse_error::SparseError;
use sparse_mat::sparse_matrix::SparseMatrix;

fn main() -> Result<(), SparseError> {
    // Assemble from (row, col, value) triplets. from_triplets has COO semantics: entries
    // that share a position are summed, which is how stamped/assembled matrices build up.
    let mut built = SparseMatrix::from_triplets(
        3,
        3,
        vec![
            (0, 0, 2.0), (0, 1, -1.0),
            (1, 0, -1.0), (1, 1, 1.0), (1, 1, 1.0), (1, 2, -1.0),
            (2, 1, -1.0), (2, 2, 2.0),
        ],
    )?;
    assert_eq!(built[(1, 1)], 2.0); // the two (1, 1) entries were summed
    built.insert(2, 2, 3.0)?; // edits are cheap while building

    // Convert once for numerical work. Products and solves return Result; the operators panic
    // on shape mismatches instead.
    let a: CsrMatrix = built.to_csr();
    let b = [1.0, 0.0, 2.0];
    let x = a.to_csc().solve(&b)?; // LU: AMD-ordered, threshold partial pivoting
    let residual = a.mul_vec(&x)?;
    assert!(residual.iter().zip(&b).all(|(r, b)| (r - b).abs() < 1e-12));

    // a is symmetric positive definite: factor it as P A Pᵀ = L Lᵀ, ordered to limit fill-in
    let chol = a.cholesky()?;
    let l = chol.l();
    let lt: CsrMatrix = l.clone().transpose().to_csr(); // transpose is O(1), as a CscMatrix
    let reordered = a.permute_symmetric(chol.permutation())?;
    assert!((&(l * &lt) - &reordered).norm_frobenius() < 1e-12);
    assert_eq!(chol.solve(&b)?.len(), 3);

    // The pattern is fixed: values change in place, and a zero stays stored
    let mut refilled = a.clone();
    refilled.values_mut()[0] = 0.0;
    assert_eq!(refilled.num_stored(), a.num_stored());

    let shifted = &a + &CsrMatrix::identity(3);
    println!("{} stored entries, trace {}", (&a * &a).num_stored(), shifted.trace()?);
    Ok(())
}
```

## What it provides

**`SparseMatrix`: building and editing**
- `new`, `empty_with_shape`, `identity`, `from_diagonal`.
- `from_triplets`, which sums duplicates (COO semantics), `from_dense`, `From<[[T; C]; R]>` and
  `TryFrom<Vec<Vec<T>>>`.
- Element edits with `insert`, `insert_triplets` (duplicates overwrite), `clear_at` and
  `prune_zeroes`; checked reads with `peek_at`, and `matrix[(row, col)]` indexing, which panics
  out of bounds.
- A zero is never stored, so `num_nonzero` counts real entries.
- `iter()` over stored entries in row-major order, `to_dense`, `create_transpose` and
  `transpose_inplace`.
- Matrix Market coordinate files through `from_matrix_market`, `read_matrix_market`,
  `write_matrix_market` and `to_matrix_market`. Reading supports `real`, `integer` and `pattern`
  fields with `general`, `symmetric` or `skew-symmetric` symmetry.
- `to_csr()` and `to_csc()` for numerical work.

**`CompressedMatrix` (`CsrMatrix`, `CscMatrix`): numerical work**
- Built with `from_triplets` (every given position is stored, even if its values cancel),
  `from_dense`, `identity`, `from_diagonal`, `from_sparse`, or `try_from_raw_parts`, which
  validates raw arrays. `into_raw_parts`, `pointers`, `indices` and `values` expose the arrays.
- A fixed pattern: `values_mut` changes values in place, and a stored entry stays stored when
  its value is zero. `num_stored` counts the pattern and `num_nonzero` the nonzero values.
- Slots in `values`, for refilling without searching: `outer(i)` and `outer_mut(i)` borrow a
  CSR row or CSC column's indices and values without copying, `position(row, col)` finds an
  entry's slot by binary search, and `diagonal_positions` gives the diagonal's slots, computed
  once and cached, since the pattern can't change.
- **Assembly maps:** `from_triplets_mapped` builds the same matrix as `from_triplets` along
  with an `AssemblyMap` of each triplet's slot, and `refill` reassembles new values for the
  same triplets in O(triplets), with no sorting or searching, bit for bit as `from_triplets`
  would. For matrices rebuilt with new values, such as a Jacobian at every Newton iteration.
- `transpose` returns the same arrays in the other orientation, in O(1); `to_csr` and `to_csc`
  reorient explicitly, in O(nnz); `to_sparse` converts back, dropping stored zeros.
- Equality is mathematical: the same shape and nonzero entries, whatever the orientation or
  stored zeros, also against a `SparseMatrix`.
- **Access:** `row(i)` and `row_iter` on CSR, `col(j)` and `col_iter` on CSC (copies of the
  slices `outer` borrows), and `iter`, `diagonal`, `submatrix`, `select` (arbitrary index lists in any order),
  `permute` and `permute_symmetric` (P A Pᵀ) in either orientation.
- **Products:** `mul_vec` (`&a * &x`), `vec_mul` (`&x * &a`), `mul_dense` and sparse × sparse
  `matmul` (`&a * &b`, Gustavson's algorithm, keeping every reachable position). Sums walk
  sorted storage, so results are deterministic, and CSR and CSC give bit-for-bit equal results.
- **Arithmetic:** `+`, `-` (storing the union of the patterns), unary `-`, `+=`, `-=`, `*=` and
  scalar `*`, `try_add`, `try_sub`, `scale`, `hadamard`, `map`/`apply` over stored entries, and
  `scale_rows`/`scale_cols` for diag(v)·A and A·diag(v). `kron`, `hstack` and `vstack`.
- **Complex matrices:** `conj`, `adjoint` (in the other orientation, like `transpose`),
  `real`/`imag`, and `from_parts`, generic over `ComplexFloat`.
- **Properties:** `is_square`, `is_symmetric`, `density`, `sum`, `trace`, and `norm_frobenius`,
  `norm_1` and `norm_inf` (float element types).
- **Solvers** (float element types): `solve_lower_triangular` and `solve_upper_triangular` in
  either orientation, and `conjugate_gradient`; on `CsrMatrix`, `cholesky()` (a `Cholesky` of
  P A Pᵀ, ordered by AMD, with `solve`, `l` and `permutation`) and `solve_spd`.
- **LU** (real or complex element types): on `CscMatrix`, `lu()` and `solve`, Gilbert and
  Peierls's left-looking LU of P A Q, AMD-ordered, with threshold partial pivoting (0.1 by
  default; `LuOptions`). The `Lu` has `solve`, `determinant`, `l`, `u`, `permutation` and
  `column_permutation`. `analyze_lu(ordering)` separates the ordering from `factor`.
- **Orderings and symbolic analysis:** `ordering(Ordering::Amd)` (approximate minimum degree)
  and `ordering(Ordering::Rcm)` (reverse Cuthill–McKee), on the pattern of A + Aᵀ in either
  orientation, with `bandwidth` and `profile` to measure the latter;
  `analyze_cholesky(ordering)` works out the elimination tree and the factor's exact pattern
  once, so `factor` can be repeated for every matrix with that pattern.

**Standard traits:** `Clone`, `Debug`, `PartialEq`, `Display` and `Default` on both types;
`IntoIterator` for `SparseMatrix`.

## Element types

Both types accept any `Scalar`, meaning any `Copy + PartialEq + num_traits::Num` type. That
covers `f32`, `f64`, every integer type and `num_complex::Complex` (re-exported as
`sparse_mat::Complex`). Arithmetic, products and conversions work for all of them. LU works
for any `num_complex::ComplexFloat` type: `f32`, `f64` and their complex numbers. The norms and
the other solvers need `num_traits::Float`, so they cover `f32` and `f64` but not complex
numbers yet; Phase 4 of the roadmap migrates them.

Rust does not apply a type parameter's default in expressions, so a call with nothing to infer
`T` from needs it spelled out: `SparseMatrix::<f64>::new()`.

## Storage model and performance

- **`SparseMatrix`** keeps a hash map keyed by `(row, col)`, so single-element edits are O(1).
  `iter()` sorts the entries when called. It works for any shape, even dimensions near
  `usize::MAX`.
- **`CompressedMatrix`** keeps three arrays: one pointer per row (CSR) or column (CSC) plus one,
  and an index and a value per stored entry. A pointer per row or column costs memory however
  sparse the matrix is, so a dimension too large for that panics on conversion with a message
  saying so.
- **Conversions are explicit.** `to_csr`, `to_csc` and reorienting cost O(nnz) plus sorting;
  `transpose` is free. Nothing is rebuilt behind your back, so convert once and keep the
  compressed matrix for repeated work.

Known performance limits:
- `to_dense` and `mul_dense` allocate `rows × cols`.

## Errors

Fallible operations return `Result<_, SparseError>`. Its variants name the problem and carry
the offending values, for example `OutOfBounds { row, col, shape }`,
`ShapeMismatch { left, right }`, `SingularMatrix { pivot }` and `Parse { line, message }`. The
enum is `#[non_exhaustive]`.

Operators (`+`, `*`, `[]`, …) cannot return `Result`, so they panic on invalid input. Each has a
checked method alongside it (`try_add`, `matmul`, `mul_vec`, `peek_at`, …).

## Optional features

| Feature | Adds |
|---|---|
| `serde` | `Serialize`/`Deserialize` for `SparseMatrix` as `{"shape": [r, c], "entries": [[row, col, value], ...]}`, with entries in row-major order and validated like `from_triplets`; and for `CompressedMatrix` as its arrays with the orientation, validated like `try_from_raw_parts` |
| `rayon` | `par_mul_vec` on `CsrMatrix`, `par_vec_mul` on `CscMatrix` and `par_matmul` on both, which return results bit-for-bit equal to the sequential versions |

## Development

- **MSRV:** Rust 1.87, set for the whole workspace in the root `Cargo.toml`.
- **Tests:** in `tests/`. The products, factorizations and Kronecker products are checked
  against dense reference implementations and algebraic identities. The Matrix Market reader
  is fuzzed with truncated, malformed and oversized input.
- **API documentation:** `cargo doc -p sparse_mat --open --all-features`.

See the [workspace README](../README.md) for `./check.sh`, CI and the rest of the docs.

## License

Licensed under either of the Apache License, Version 2.0 ([`LICENSE-APACHE`](LICENSE-APACHE))
or the MIT license ([`LICENSE-MIT`](LICENSE-MIT)), at your option. See
[`AI-DISCLAIMER.md`](AI-DISCLAIMER.md) for how this licensing relates to the use of AI in
developing the crate.

## References

The published algorithms behind this crate, and the license terms of related software, are
listed in the [workspace README](../README.md#references).
