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
