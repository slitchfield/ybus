// The GO Challenge 1 case.raw reader on the hand-written fixture: every field type read, the
//   line each record came from, and the text variants G2 found in the data (CRLF, mixed line
//   endings, Q in place of the last terminators, no Q at all) reading the same.

use std::fs;
use std::path::PathBuf;

use ybus_powerflow::case::go::RawCase;
use ybus_powerflow::case::BusType;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/go/small/case.raw")
}

/// The fixture's text, committed with LF line endings
fn fixture() -> String {
    fs::read_to_string(fixture_path()).expect("the fixture is committed")
}

/// The fixture's lines, without their endings
fn fixture_lines() -> Vec<String> {
    fixture().lines().map(str::to_string).collect()
}

#[test]
fn the_fixture_reads_field_by_field() {
    let case = RawCase::parse(&fixture()).unwrap();
    assert_eq!(case.sbase_mva, 100.0);

    // Buses: names trimmed, IDE as a bus type, both voltage limit pairs
    assert_eq!(case.buses.len(), 4);
    let bus = &case.buses[0];
    assert_eq!((bus.source_line, bus.bus_label), (4, 1));
    assert_eq!((bus.name.as_str(), bus.base_kv), ("BUS ONE", 230.0));
    assert_eq!((bus.bus_type, bus.area, bus.zone), (BusType::Ref, 1, 1));
    assert_eq!((bus.vm_pu, bus.va_deg), (1.04, 0.0));
    let bus = &case.buses[2];
    assert_eq!(
        (bus.name.as_str(), bus.bus_type, bus.zone),
        ("", BusType::PQ, 2)
    );
    assert_eq!(
        (bus.nvhi_pu, bus.nvlo_pu, bus.evhi_pu, bus.evlo_pu),
        (1.1, 0.9, 1.15, 0.85)
    );
    assert_eq!(case.buses[1].bus_type, BusType::PV);
    assert_eq!(case.buses[1].name, "");

    // Two loads and two fixed shunts on bus 3, with padded IDs trimmed
    let ids: Vec<_> = case
        .loads
        .iter()
        .map(|l| (l.bus_label, l.id.as_str(), l.status))
        .collect();
    assert_eq!(ids, [(3, "1", true), (3, "2", true), (4, "L", false)]);
    assert_eq!((case.loads[1].pl_mw, case.loads[1].ql_mvar), (10.5, -3.25));
    assert_eq!(case.loads[0].source_line, 9);
    let shunts: Vec<_> = case
        .fixed_shunts
        .iter()
        .map(|s| (s.id.as_str(), s.status, s.gl_mw, s.bl_mvar))
        .collect();
    assert_eq!(shunts, [("1", true, 0.0, 15.0), ("2", false, 1.5, -5.0)]);

    // Generators, one out of service and one with MBASE 0, in exponent notation in places
    assert_eq!(case.generators.len(), 3);
    let g = &case.generators[0];
    assert_eq!(
        (g.bus_label, g.id.as_str(), g.pg_mw, g.qg_mvar),
        (1, "1", 80.0, 10.0)
    );
    assert_eq!((g.qt_mvar, g.qb_mvar, g.mbase_mva), (100.0, -50.0, 100.0));
    assert_eq!((g.status, g.pt_mw, g.pb_mw), (true, 200.0, 10.0));
    assert_eq!(case.generators[1].mbase_mva, 0.0);
    let g = &case.generators[2];
    assert_eq!((g.id.as_str(), g.status, g.source_line), ("99", false, 18));

    // Lines: a letter circuit, one out of service with no resistance
    let l = &case.lines[0];
    assert_eq!((l.from_label, l.to_label, l.circuit.as_str()), (1, 2, "1"));
    assert_eq!((l.r_pu, l.x_pu, l.b_pu), (0.01, 0.1, 0.02));
    assert_eq!(
        (l.rate_a_mva, l.rate_b_mva, l.rate_c_mva),
        (250.0, 0.0, 300.0)
    );
    assert_eq!((l.gi_pu, l.bi_pu, l.gj_pu, l.bj_pu), (0.0, 0.0, 0.0, 0.0));
    let l = &case.lines[1];
    assert_eq!((l.circuit.as_str(), l.status, l.r_pu), ("TL", false, 0.0));

    // Transformers: four lines each, a magnetizing susceptance, an off-nominal ratio and a
    //   shift, and one with CZ 2 out of service
    assert_eq!(case.transformers.len(), 2);
    let t = &case.transformers[0];
    assert_eq!((t.source_line, t.from_label, t.to_label), (24, 2, 3));
    assert_eq!((t.circuit.as_str(), t.cw, t.cz, t.cm), ("1", 1, 1, 1));
    assert_eq!((t.mag1_pu, t.mag2_pu, t.status), (0.0, -0.0025, true));
    assert_eq!((t.r12_pu, t.x12_pu, t.sbase12_mva), (0.002, 0.08, 100.0));
    assert_eq!((t.windv1, t.nomv1_kv, t.ang1_deg), (1.025, 230.0, -3.0));
    assert_eq!((t.rata1_mva, t.ratb1_mva, t.ratc1_mva), (200.0, 0.0, 240.0));
    assert_eq!((t.windv2, t.nomv2_kv), (1.0, 115.0));
    let t = &case.transformers[1];
    assert_eq!(
        (t.source_line, t.circuit.as_str(), t.cz, t.status),
        (28, "2", 2, false)
    );

    // The switched shunt, after the unused sections, which hold records
    assert_eq!(case.switched_shunts.len(), 1);
    let s = &case.switched_shunts[0];
    assert_eq!(
        (s.source_line, s.bus_label, s.status, s.binit_mvar),
        (48, 4, true, 12.0)
    );
    let mut blocks = [(0, 0.0); 8];
    blocks[0] = (2, 10.0);
    blocks[1] = (1, -5.0);
    assert_eq!(s.blocks, blocks);
}

#[test]
fn line_endings_and_file_endings_read_the_same() {
    let want = RawCase::parse(&fixture()).unwrap();
    let lines = fixture_lines();
    let join = |lines: &[String], ending: &str| lines.join(ending) + ending;

    // CRLF throughout, and LF on lines 2 and 3 only, as 60 of the data's files have them
    let crlf = join(&lines, "\r\n");
    let mut mixed = String::new();
    for (i, line) in lines.iter().enumerate() {
        mixed += line;
        mixed += if i == 1 || i == 2 { "\n" } else { "\r\n" };
    }
    // Q in place of the GNE and induction machine terminators, as 60 files end
    let q_early: Vec<String> = lines
        .iter()
        .filter(|l| !(l.contains("END OF GNE DEVICE DATA") || l.contains("END OF INDUCTION")))
        .cloned()
        .collect();
    assert_eq!(q_early.len(), lines.len() - 2);
    // No Q, and the end of the file right after the switched shunt section
    let no_q = &lines[..lines.len() - 1];
    let shunts_end = lines
        .iter()
        .position(|l| l.starts_with("0 / END OF SWITCHED SHUNT DATA"))
        .unwrap();
    let ends_early = &lines[..=shunts_end];

    for (what, text) in [
        ("CRLF", crlf),
        ("mixed endings", mixed),
        ("Q early", join(&q_early, "\n")),
        ("no Q", join(no_q, "\n")),
        ("ends after the switched shunts", join(ends_early, "\n")),
        ("no final newline", lines.join("\n")),
    ] {
        assert_eq!(RawCase::parse(&text), Ok(want.clone()), "{what}");
    }
}

#[test]
fn path_and_bytes_read_as_the_text_does() {
    let want = RawCase::parse(&fixture()).unwrap();
    assert_eq!(RawCase::from_path(fixture_path()), Ok(want.clone()));
    assert_eq!(RawCase::from_bytes(fixture().as_bytes()), Ok(want.clone()));
    // A byte order mark is ignored, as for MATPOWER files
    let mut bom = b"\xEF\xBB\xBF".to_vec();
    bom.extend_from_slice(fixture().as_bytes());
    assert_eq!(RawCase::from_bytes(&bom), Ok(want));
    let err = RawCase::from_path("no/such/case.raw").unwrap_err();
    assert!(err
        .to_string()
        .starts_with("can't read case file no/such/case.raw:"));
}
