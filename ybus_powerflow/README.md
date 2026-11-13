# ybus_powerflow

Power flow for electrical grids, validated against
[MATPOWER](https://matpower.org)'s own results:

- **Cases:** a model of a power system (buses, generators, branches and generator costs), read
  from MATPOWER case files (`Case::from_path`, `Case::parse`).
- **DC power flow** (`dc::dc_power_flow`), as MATPOWER's `rundcpf` computes it, with studies
  that reuse one factorization: many injections at once, single-branch outages, and PTDF and
  LODF sensitivities (`dc::DcStudy`, `dc::DcSensitivity`).
- **AC power flow** (`ac::ac_power_flow`), as MATPOWER's `runpf` computes it: by Newton's
  method or the fast-decoupled XB and BX methods, optionally enforcing generators' reactive
  limits.
- **DC and AC optimal power flow** (`opf::dc::dc_opf`, `opf::ac::ac_opf`), as MATPOWER's
  `rundcopf` and `runopf` compute them with MIPS, or by a predictor-corrector method. The
  solvers are [`ybus_ipm`](https://crates.io/crates/ybus_ipm)'s.
- **Islands and bridges** of a network (`islands::find_islands`, `islands::find_bridges`).

The linear algebra is [`ybus_sparse`](https://crates.io/crates/ybus_sparse)'s.

## Status

0.1, the first release. Each solver is tested against results computed by MATPOWER, on its
own case files, but the API is young: expect breaking changes in the 0.x releases, each one a minor
version bump.

## Quick start

```toml
[dependencies]
ybus_powerflow = "0.1"
```

A three-bus case, written as a MATPOWER case file would be, solved by the DC and the AC power
flow:

```rust
use ybus_powerflow::ac::ac_power_flow;
use ybus_powerflow::dc::dc_power_flow;
use ybus_powerflow::{Case, PowerflowError};

const CASE: &str = "
function mpc = three_bus
mpc.version = '2';
mpc.baseMVA = 100;
%% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
    1 3   0  0 0 0 1 1 0 135 1 1.1 0.9;
    2 2  50 10 0 0 1 1 0 135 1 1.1 0.9;
    3 1 100 30 0 0 1 1 0 135 1 1.1 0.9;
];
%% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin, and eleven unused columns
mpc.gen = [
    1  0 0 100 -100 1    100 1 300 0 0 0 0 0 0 0 0 0 0 0 0;
    2 80 0  50  -50 1.01 100 1 200 0 0 0 0 0 0 0 0 0 0 0 0;
];
%% fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
    1 2 0.02 0.06 0.03 0 0 0 0 0 1 -360 360;
    2 3 0.03 0.09 0.02 0 0 0 0 0 1 -360 360;
    1 3 0.01 0.05 0.02 0 0 0 0 0 1 -360 360;
];
";

fn main() -> Result<(), PowerflowError> {
    let case = Case::parse(CASE)?; // or Case::from_path("case9.m")

    // DC: lossless, so the reference generator makes up the 150 MW of load less bus 2's 80
    let dc = dc_power_flow(&case)?;
    assert!((dc.pg_mw[0] - 70.0).abs() < 1e-9);

    // AC: the reference generator also supplies the losses
    let ac = ac_power_flow(&case)?;
    assert!(ac.pg_mw[0] > dc.pg_mw[0]);
    println!(
        "{} Newton iterations, |V| at bus 3: {:.4} p.u.",
        ac.iterations, ac.vm_pu[2]
    );
    Ok(())
}
```

## Features

| Feature | Adds |
|---|---|
| `rayon` | Parallel batches of outages and of PTDF and LODF columns (`DcStudy::par_solve_outages`, `DcSensitivity::par_ptdf_columns` and `par_lodf_columns`), and `ybus_sparse`'s `rayon` feature |
| `cuda` | **Experimental.** The DC studies on a GPU, through [`ybus_cuda`](https://crates.io/crates/ybus_cuda): batched solves, single-branch outages and PTDF columns (`DcStudy::gpu`, `DcSensitivity::gpu`). It builds without the CUDA toolkit; running it needs a CUDA 13 driver and cuSPARSE. |

## Development

- **MSRV:** Rust 1.88, set for the whole workspace in the root `Cargo.toml`.
- **Tests:** in the repository's `ybus_powerflow/tests/`. They compare every power flow and
  optimal power flow with results computed by MATPOWER, on its own case files, which
  `scripts/fetch-matpower.sh` downloads. Neither the tests nor those results are in the
  published package.
- **API documentation:** `cargo doc -p ybus_powerflow --open --all-features`.

See the [workspace README](https://github.com/slitchfield/ybus#readme) for `./check.sh`, CI,
the benchmarks and the rest of the docs.

## License

Licensed under either of the Apache License, Version 2.0
([`LICENSE-APACHE`](https://github.com/slitchfield/ybus/blob/main/LICENSE-APACHE)) or the MIT
license ([`LICENSE-MIT`](https://github.com/slitchfield/ybus/blob/main/LICENSE-MIT)), at your
option. See
[`AI-DISCLAIMER.md`](https://github.com/slitchfield/ybus/blob/main/AI-DISCLAIMER.md) for how
this licensing relates to the use of AI in developing the crate.

The power flows reimplement MATPOWER's behaviour in Rust, written from its documentation and
with reference to its BSD-licensed source, rather than translating it. The workspace README's
[references](https://github.com/slitchfield/ybus#references) name the MATPOWER functions each
part follows, and the terms of MATPOWER's case files, which the repository doesn't include.
