//! Product and sum plans against `matmul` and `try_add`: the same pattern, and every value bit
//! for bit, on random real and complex matrices in both orientations. The values include
//! stored zeros, `-0.0`, cancellations and magnitudes far apart, so a term added in another
//! order, or a zero added where a value should be copied, changes some bits.

use ybus_sparse::compressed_matrix::{CscMatrix, CsrMatrix, Orientation};
use ybus_sparse::sparse_error::SparseError;
use ybus_sparse::{Complex, CompressedMatrix};

type C64 = Complex<f64>;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }

    /// A value from a mix that makes the order of additions matter: zeros of both signs, ±1,
    /// and magnitudes from 1e-8 to 1e16
    fn value(&mut self) -> f64 {
        match self.next() % 8 {
            0 => 0.0,
            1 => -0.0,
            2 => 1.0,
            3 => -1.0,
            4 => 1e16,
            5 => -1e16,
            6 => 1e-8 * (self.next() % 1000) as f64,
            _ => (self.next() % 2001) as f64 / 1000.0 - 1.0,
        }
    }
}

/// Triplets of a random `rows × cols` pattern with about `percent` of the positions stored,
/// with values from `value`. Every given position is stored, zeros included.
fn random_triplets<T>(
    rows: usize,
    cols: usize,
    percent: u64,
    rng: &mut Lcg,
    mut value: impl FnMut(&mut Lcg) -> T,
) -> Vec<(usize, usize, T)> {
    let mut triplets = vec![];
    for i in 0..rows {
        for j in 0..cols {
            if rng.next() % 100 < percent {
                triplets.push((i, j, value(rng)));
            }
        }
    }
    triplets
}

fn real(rng: &mut Lcg) -> f64 {
    rng.value()
}

fn complex(rng: &mut Lcg) -> C64 {
    C64::new(rng.value(), rng.value())
}

trait Bits {
    fn bits(&self) -> Vec<u64>;
}

impl Bits for [f64] {
    fn bits(&self) -> Vec<u64> {
        self.iter().map(|v| v.to_bits()).collect()
    }
}

impl Bits for [C64] {
    fn bits(&self) -> Vec<u64> {
        self.iter()
            .flat_map(|v| [v.re.to_bits(), v.im.to_bits()])
            .collect()
    }
}

/// The same pattern, and the same values bit for bit
fn assert_identical<T, O: Orientation>(
    planned: &CompressedMatrix<T, O>,
    direct: &CompressedMatrix<T, O>,
) where
    T: ybus_sparse::Scalar,
    [T]: Bits,
{
    assert_eq!(planned.shape(), direct.shape());
    assert_eq!(planned.pointers(), direct.pointers());
    assert_eq!(planned.indices(), direct.indices());
    assert_eq!(planned.values().bits(), direct.values().bits());
}

/// The shapes `(rows, inner, cols)` and fill percentages tried: square and rectangular, with
/// empty rows and columns in the sparse ones
const SHAPES: [(usize, usize, usize, u64); 6] = [
    (1, 1, 1, 100),
    (3, 5, 2, 60),
    (8, 8, 8, 30),
    (12, 4, 20, 25),
    (30, 30, 30, 8),
    (40, 25, 15, 3),
];

fn check_products<T, O: Orientation>(value: impl Fn(&mut Lcg) -> T + Copy)
where
    T: ybus_sparse::Scalar,
    [T]: Bits,
{
    for seed in 0..20u64 {
        for (k, &(rows, inner, cols, fill)) in SHAPES.iter().enumerate() {
            let mut rng = Lcg(seed * 31 + k as u64);
            let a: CompressedMatrix<T, O> = CompressedMatrix::from_triplets(
                rows,
                inner,
                random_triplets(rows, inner, fill, &mut rng, value),
            )
            .unwrap();
            let b: CompressedMatrix<T, O> = CompressedMatrix::from_triplets(
                inner,
                cols,
                random_triplets(inner, cols, fill, &mut rng, value),
            )
            .unwrap();
            let plan = a.plan_product(&b).unwrap();
            let direct = a.matmul(&b).unwrap();
            assert_identical(&plan.multiply(&a, &b).unwrap(), &direct);
            assert_eq!(plan.num_stored(), direct.num_stored());

            // The same plan with new values for both operands
            let (mut a2, mut b2) = (a.clone(), b.clone());
            a2.values_mut()
                .iter_mut()
                .for_each(|v| *v = value(&mut rng));
            b2.values_mut()
                .iter_mut()
                .for_each(|v| *v = value(&mut rng));
            assert_identical(&plan.multiply(&a2, &b2).unwrap(), &a2.matmul(&b2).unwrap());

            // Into a reused buffer that holds the previous values
            let mut out = plan.multiply(&a, &b).unwrap().values().to_vec();
            plan.multiply_values(a2.values(), b2.values(), &mut out)
                .unwrap();
            assert_eq!(out.bits(), a2.matmul(&b2).unwrap().values().bits());
        }
    }
}

fn check_sums<T, O: Orientation>(value: impl Fn(&mut Lcg) -> T + Copy)
where
    T: ybus_sparse::Scalar,
    [T]: Bits,
{
    for seed in 0..20u64 {
        for (k, &(rows, _, cols, fill)) in SHAPES.iter().enumerate() {
            let mut rng = Lcg(seed * 37 + k as u64);
            let a: CompressedMatrix<T, O> = CompressedMatrix::from_triplets(
                rows,
                cols,
                random_triplets(rows, cols, fill, &mut rng, value),
            )
            .unwrap();
            let b: CompressedMatrix<T, O> = CompressedMatrix::from_triplets(
                rows,
                cols,
                random_triplets(rows, cols, fill, &mut rng, value),
            )
            .unwrap();
            let plan = a.plan_sum(&b).unwrap();
            assert_identical(&plan.add(&a, &b).unwrap(), &a.try_add(&b).unwrap());

            let (mut a2, mut b2) = (a.clone(), b.clone());
            a2.values_mut()
                .iter_mut()
                .for_each(|v| *v = value(&mut rng));
            b2.values_mut()
                .iter_mut()
                .for_each(|v| *v = value(&mut rng));
            assert_identical(&plan.add(&a2, &b2).unwrap(), &a2.try_add(&b2).unwrap());
            // b + a's plan is a different plan: the operands' roles swap
            assert_identical(
                &b2.plan_sum(&a2).unwrap().add(&b2, &a2).unwrap(),
                &b2.try_add(&a2).unwrap(),
            );
        }
    }
}

#[test]
fn real_products_match_matmul_bit_for_bit() {
    check_products::<f64, ybus_sparse::compressed_matrix::Csr>(real);
    check_products::<f64, ybus_sparse::compressed_matrix::Csc>(real);
}

#[test]
fn complex_products_match_matmul_bit_for_bit() {
    check_products::<C64, ybus_sparse::compressed_matrix::Csr>(complex);
    check_products::<C64, ybus_sparse::compressed_matrix::Csc>(complex);
}

#[test]
fn real_sums_match_try_add_bit_for_bit() {
    check_sums::<f64, ybus_sparse::compressed_matrix::Csr>(real);
    check_sums::<f64, ybus_sparse::compressed_matrix::Csc>(real);
}

#[test]
fn complex_sums_match_try_add_bit_for_bit() {
    check_sums::<C64, ybus_sparse::compressed_matrix::Csr>(complex);
    check_sums::<C64, ybus_sparse::compressed_matrix::Csc>(complex);
}

#[test]
fn the_order_of_terms_is_matmuls() {
    // (1e16 + 1) - 1e16 is 0 in this order and 1 in others, so only matmul's order passes
    let a = CsrMatrix::from_triplets(1, 3, vec![(0, 0, 1e16), (0, 1, 1.0), (0, 2, -1e16)]).unwrap();
    let b = CsrMatrix::from_triplets(3, 1, vec![(0, 0, 1.0), (1, 0, 1.0), (2, 0, 1.0)]).unwrap();
    let product = a.plan_product(&b).unwrap().multiply(&a, &b).unwrap();
    assert_eq!(product.values(), a.matmul(&b).unwrap().values());
    assert_eq!(product.values(), &[0.0]);
}

#[test]
fn a_copied_negative_zero_stays_negative() {
    // try_add copies an entry only one operand stores, so -0.0 stays -0.0; adding 0.0 would
    //   make it +0.0
    let a = CsrMatrix::from_triplets(1, 2, vec![(0, 0, -0.0)]).unwrap();
    let b = CsrMatrix::from_triplets(1, 2, vec![(0, 1, -0.0)]).unwrap();
    let sum = a.plan_sum(&b).unwrap().add(&a, &b).unwrap();
    assert_eq!(sum.values().bits(), [(-0.0f64).to_bits(); 2]);
    // And a product's entry starts at +0.0, as matmul's does: 0 + (-0.0 · 1) is +0.0
    let one = CsrMatrix::from_triplets(2, 1, vec![(0, 0, 1.0)]).unwrap();
    let product = a.plan_product(&one).unwrap().multiply(&a, &one).unwrap();
    assert_eq!(
        product.values().bits(),
        a.matmul(&one).unwrap().values().bits()
    );
    assert_eq!(product.values().bits(), [0.0f64.to_bits()]);
}

#[test]
fn reachable_positions_are_stored_even_where_terms_cancel() {
    let a = CsrMatrix::from_triplets(1, 2, vec![(0, 0, 1.0), (0, 1, -1.0)]).unwrap();
    let b = CsrMatrix::from_triplets(2, 1, vec![(0, 0, 2.0), (1, 0, 2.0)]).unwrap();
    let product = a.plan_product(&b).unwrap().multiply(&a, &b).unwrap();
    assert_eq!(product.num_stored(), 1);
    assert_eq!(product.values(), &[0.0]);
}

#[test]
fn mismatched_operands_are_errors() {
    let a = CsrMatrix::from_triplets(2, 3, vec![(0, 0, 1.0), (1, 2, 2.0)]).unwrap();
    let b = CsrMatrix::from_triplets(3, 2, vec![(0, 1, 1.0), (2, 0, 3.0)]).unwrap();
    assert!(matches!(
        a.plan_product(&a),
        Err(SparseError::ShapeMismatch { .. })
    ));
    assert!(matches!(
        a.plan_sum(&b),
        Err(SparseError::ShapeMismatch { .. })
    ));

    let product = a.plan_product(&b).unwrap();
    // Another pattern of the same shape
    let other = CsrMatrix::from_triplets(2, 3, vec![(0, 0, 1.0)]).unwrap();
    assert_eq!(
        product.multiply(&other, &b),
        Err(SparseError::PatternMismatch)
    );
    let mut out = vec![7.0; product.num_stored()];
    assert_eq!(
        product.multiply_values(&[1.0], b.values(), &mut out),
        Err(SparseError::LengthMismatch {
            expected: 2,
            found: 1
        })
    );
    assert_eq!(
        out,
        vec![7.0; product.num_stored()],
        "unchanged on an error"
    );
    assert!(matches!(
        product.multiply_values(a.values(), b.values(), &mut [0.0; 9]),
        Err(SparseError::LengthMismatch { .. })
    ));

    let sum = a.plan_sum(&a).unwrap();
    assert_eq!(sum.add(&other, &a), Err(SparseError::PatternMismatch));
    let c: CscMatrix<f64> = a.to_csc();
    assert!(c.plan_sum(&c).is_ok());
}
