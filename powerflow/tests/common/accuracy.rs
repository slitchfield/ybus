//! How accurately a DC power flow's angles solve the system they came from, a measure that
//! needs no reference results and doesn't depend on the case's size or conditioning.

use std::collections::HashMap;
use std::f64::consts::PI;

use powerflow::case::Case;
use powerflow::dc::DcModel;

/// The normwise backward error of `va_deg` (angles in degrees, one per case bus row) as a
/// solution of the DC system `model` builds for `case`:
///
/// η = ‖P′ − B′θ‖∞ / (‖B′‖∞ ‖θ‖∞ + ‖P′‖∞)
///
/// over the PV and PQ buses, with the reference angles taken from `va_deg`. It is the smallest
/// relative change to B′ and P′ that makes θ an exact solution, so a backward-stable solve
/// keeps it within a small multiple of `f64::EPSILON`, however ill-conditioned B′ is. P′ is
/// derived here, independently of `DcModel::solve`: (Pg − Pd − Gs) / baseMVA over in-service
/// generators, less the phase shifters' injections and the reference buses' coupling.
///
/// 0.0 when every connected bus is a reference bus. NaN if an angle is NaN.
pub fn backward_error(case: &Case, model: &DcModel, va_deg: &[f64]) -> f64 {
    let rows = model.bus_rows();
    let unknown: Vec<usize> = model
        .pv_buses()
        .iter()
        .chain(model.pq_buses())
        .copied()
        .collect();
    if unknown.is_empty() {
        return 0.0;
    }

    let mut internal = vec![usize::MAX; case.bus.len()];
    for (index, &row) in rows.iter().enumerate() {
        internal[row] = index;
    }
    let by_label: HashMap<usize, usize> = case.bus_indices();
    let base_mva = case.base_mva;
    let mut pbus: Vec<f64> = rows
        .iter()
        .zip(model.pbusinj_pu())
        .map(|(&row, pbusinj)| -(case.bus[row].pd_mw + case.bus[row].gs_mw) / base_mva - pbusinj)
        .collect();
    for &row in model.gen_rows() {
        let generator = &case.generator[row];
        pbus[internal[by_label[&generator.bus_label]]] += generator.pg_mw / base_mva;
    }

    let radians = |index: usize| va_deg[rows[index]] * PI / 180.0;
    let theta: Vec<f64> = unknown.iter().map(|&index| radians(index)).collect();
    let theta_ref: Vec<f64> = model
        .ref_buses()
        .iter()
        .map(|&index| radians(index))
        .collect();
    let bprime = model
        .bbus()
        .select(&unknown, &unknown)
        .expect("bus indices are in range");
    let coupling = model
        .bbus()
        .select(&unknown, model.ref_buses())
        .expect("bus indices are in range")
        .mul_vec(&theta_ref)
        .expect("one angle per reference bus");
    let p: Vec<f64> = unknown
        .iter()
        .zip(&coupling)
        .map(|(&index, coupling)| pbus[index] - coupling)
        .collect();
    let b_theta = bprime.mul_vec(&theta).expect("one angle per bus");
    let residual: Vec<f64> = p.iter().zip(&b_theta).map(|(p, bt)| p - bt).collect();

    let mut row_sums = vec![0.0; unknown.len()];
    for (row, _, value) in bprime.iter() {
        row_sums[row] += value.abs();
    }
    let scale = largest_magnitude(&row_sums) * largest_magnitude(&theta) + largest_magnitude(&p);
    let residual = largest_magnitude(&residual);
    if residual == 0.0 {
        0.0
    } else {
        residual / scale
    }
}

/// The largest absolute value, or NaN if any value is NaN, which `f64::max` would skip.
fn largest_magnitude(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |largest, value| {
        if value.is_nan() || value.abs() > largest {
            value.abs()
        } else {
            largest
        }
    })
}
