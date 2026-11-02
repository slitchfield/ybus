// Tests for the AC power flow against MATPOWER's runpf results, committed for case14, case30
//   and case118 in tests/data/matpower-reference/.

mod common;
use common::reference::{committed_dir, read_ac};

#[test]
fn committed_ac_references_read_and_describe_their_cases() {
    for (name, buses, branches, gens, iterations) in [
        ("case14", 14, 20, 5, 2),
        ("case30", 30, 41, 6, 3),
        ("case118", 118, 186, 54, 3),
    ] {
        let reference = read_ac(&committed_dir(), name).unwrap();
        assert_eq!(reference.bus.len(), buses, "{name}");
        assert_eq!(reference.branch.len(), branches, "{name}");
        assert_eq!(reference.gen.len(), gens, "{name}");
        assert!(reference.success, "{name}");
        assert_eq!(reference.iterations, iterations, "{name}");
        // Voltages near 1 p.u., and a lossy network: at least one branch's flows don't cancel
        assert!(reference
            .bus
            .iter()
            .all(|&(_, vm, _)| (0.9..1.1).contains(&vm)));
        assert!(reference
            .branch
            .iter()
            .any(|(_, _, _, [pf, _, pt, _])| (pf + pt).abs() > 1e-6));
    }
}

/// The validation tolerance agreed for Phase 4.5, as for the DC power flow: every result
/// within 1e-9 of MATPOWER's (p.u., degrees, MW, MVAr).
const TOLERANCE: f64 = 1e-9;

#[test]
fn newton_raphson_matches_runpf() {
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue; // not fetched; the skip was reported on stderr
        };
        let case = powerflow::case::Case::from_path(&path).unwrap();
        let reference = read_ac(&committed_dir(), name).unwrap();
        let solution = powerflow::ac::ac_power_flow(&case).unwrap();
        let deviation = common::reference::ac_deviation(&case, &solution, &reference).unwrap();
        eprintln!(
            "{name}: {} iterations (MATPOWER {}), mismatch {:e}, deviation {deviation:?}",
            solution.iterations, reference.iterations, solution.mismatch_pu
        );
        assert_eq!(solution.iterations, reference.iterations, "{name}");
        assert!(
            deviation.largest() <= TOLERANCE,
            "{name}: {deviation:?} exceeds {TOLERANCE}"
        );
        assert!(solution.mismatch_pu < 1e-8);
        assert_eq!(
            solution.full_factorizations, 1,
            "{name}: later iterations refactor"
        );
    }
}

#[test]
fn fast_decoupled_matches_runpf() {
    use powerflow::ac::{AcMethod, AcModel, AcOptions};
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue; // not fetched; the skip was reported on stderr
        };
        let case = powerflow::case::Case::from_path(&path).unwrap();
        let model = AcModel::new(&case).unwrap();
        for (run, method) in [
            ("fdxb", AcMethod::FastDecoupledXb),
            ("fdbx", AcMethod::FastDecoupledBx),
        ] {
            let reference = common::reference::read_ac_run(&committed_dir(), name, run).unwrap();
            let solution = model.solve(&case, AcOptions::for_method(method)).unwrap();
            let deviation = common::reference::ac_deviation(&case, &solution, &reference).unwrap();
            eprintln!(
                "{name} {run}: {} iterations (MATPOWER {}), mismatch {:e}, largest deviation {:e}",
                solution.iterations,
                reference.iterations,
                solution.mismatch_pu,
                deviation.largest()
            );
            assert_eq!(solution.iterations, reference.iterations, "{name} {run}");
            assert!(
                deviation.largest() <= TOLERANCE,
                "{name} {run}: {deviation:?} exceeds {TOLERANCE}"
            );
            assert_eq!(solution.full_factorizations, 2, "B′ and B″, once each");
        }
    }
}

#[test]
fn reactive_limits_match_runpf() {
    use powerflow::ac::{AcModel, AcOptions, ReactiveLimits};
    for name in ["case14", "case30", "case118"] {
        let Some(path) = common::matpower_case(name) else {
            continue; // not fetched; the skip was reported on stderr
        };
        let case = powerflow::case::Case::from_path(&path).unwrap();
        let model = AcModel::new(&case).unwrap();
        for (run, mode) in [
            ("qlim1", ReactiveLimits::AllAtOnce),
            ("qlim2", ReactiveLimits::OneAtATime),
        ] {
            let reference = common::reference::read_ac_run(&committed_dir(), name, run).unwrap();
            let held = common::reference::read_limited(&committed_dir(), name, run).unwrap();
            let mut options = AcOptions::default();
            options.reactive_limits = mode;
            let solution = model.solve(&case, options).unwrap();
            let deviation = common::reference::ac_deviation(&case, &solution, &reference).unwrap();
            eprintln!(
                "{name} {run}: {} iterations (MATPOWER {}), held {:?} (MATPOWER {held:?}), largest \
                 deviation {:e}",
                solution.iterations,
                reference.iterations,
                solution.limited_generators,
                deviation.largest()
            );
            assert_eq!(solution.iterations, reference.iterations, "{name} {run}");
            let mut ours: Vec<usize> = solution.limited_generators.iter().map(|r| r + 1).collect();
            ours.sort_unstable();
            assert_eq!(ours, held, "{name} {run}");
            assert!(
                deviation.largest() <= TOLERANCE,
                "{name} {run}: {deviation:?} exceeds {TOLERANCE}"
            );
        }
    }
}
