// Set of tests for our initial data model parsing efforts

mod common;
use ybus_powerflow::case::{Branch, BusType, Case, GenCost};
use ybus_powerflow::powerflow_error::PowerflowError;

use std::path::{Path, PathBuf};

/// A small valid version 2 case: a reference bus and a load bus, one generator and one line.
/// Each table can be replaced to build the edge cases below.
struct Tiny {
    bus: &'static str,
    gen: &'static str,
    branch: &'static str,
    extra: &'static str,
}

impl Default for Tiny {
    fn default() -> Self {
        Tiny {
            bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
                  \t2\t1\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
            gen: "\t1\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
            branch: "\t1\t2\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
            extra: "",
        }
    }
}

impl Tiny {
    fn text(&self) -> String {
        format!(
            "function mpc = tiny\n\
             mpc.version = '2';\n\
             mpc.baseMVA = 100;\n\
             mpc.bus = [\n{}];\n\
             mpc.gen = [\n{}];\n\
             mpc.branch = [\n{}];\n\
             {}",
            self.bus, self.gen, self.branch, self.extra
        )
    }

    fn parse(&self) -> Result<Case, PowerflowError> {
        Case::parse(&self.text())
    }
}

/// The message of a `Parse` error, or a panic naming what came back instead.
fn parse_message(result: Result<Case, PowerflowError>) -> (usize, String) {
    match result {
        Err(PowerflowError::Parse { line, message }) => (line, message),
        other => panic!("expected a Parse error, got {other:?}"),
    }
}

#[test]
fn parse_bad_casepath() {
    match Case::from_path(PathBuf::from("trash")) {
        Err(PowerflowError::NoCaseFile { bad_path, reason }) => {
            assert_eq!(bad_path, PathBuf::from("trash"));
            // The I/O error's own message is kept, so the cause isn't lost
            assert!(!reason.is_empty());
        }
        other => panic!("expected NoCaseFile, got {other:?}"),
    }
}

#[test]
fn parse_case14() {
    let Some(path) = common::matpower_case("case14") else {
        return;
    };
    let case = Case::from_path(&path).unwrap();

    assert_eq!(case.base_mva, 100.0);
    assert_eq!(
        (case.bus.len(), case.generator.len(), case.branch.len()),
        (14, 5, 20)
    );

    let count = |kind| case.bus.iter().filter(|b| b.bus_type == kind).count();
    assert_eq!(
        (
            count(BusType::Ref),
            count(BusType::PV),
            count(BusType::PQ),
            count(BusType::Isolated)
        ),
        (1, 4, 9, 0)
    );

    let bus9 = &case.bus[8];
    assert_eq!((bus9.bus_label, bus9.pd_mw, bus9.bs_mvar), (9, 29.5, 19.0));
    assert_eq!(case.bus[0].va_deg, 0.0);

    // Row 8 is the 4-7 transformer; row 1 is a line, whose file tap of 0 means 1
    let transformer = &case.branch[7];
    assert_eq!((transformer.from_label, transformer.to_label), (4, 7));
    assert_eq!(transformer.tap_ratio, 0.978);
    assert_eq!(transformer.effective_tap(), 0.978);
    assert_eq!(case.branch[0].tap_ratio, 0.0);
    assert_eq!(case.branch[0].effective_tap(), 1.0);
    assert!(case.branch.iter().all(|b| b.status));

    let gen1 = &case.generator[0];
    assert_eq!((gen1.bus_label, gen1.pg_mw, gen1.status), (1, 232.4, true));
    assert_eq!(gen1.pmax_mw, 332.4);

    let costs = case.gencost.as_ref().unwrap();
    assert_eq!(costs.len(), 5);
    assert_eq!(
        costs[0],
        GenCost::Polynomial {
            startup: 0.0,
            shutdown: 0.0,
            coeffs: vec![0.0430292599, 20.0, 0.0],
        }
    );

    let names = case.busnames.as_ref().unwrap();
    assert_eq!(names.len(), 14);
    assert_eq!(names[0], "Bus 1     HV");
}

#[test]
fn parse_all_test_cases() {
    // The cases the power flow tests use, with their bus, generator and branch counts
    for (name, sizes) in [
        ("case14", (14, 5, 20)),
        ("case30", (30, 6, 41)),
        ("case118", (118, 54, 186)),
        ("case2383wp", (2383, 327, 2896)),
    ] {
        let Some(path) = common::matpower_case(name) else {
            return;
        };
        let case = Case::from_path(&path).unwrap_or_else(|err| panic!("{name}: {err}"));
        assert_eq!(
            (case.bus.len(), case.generator.len(), case.branch.len()),
            sizes,
            "{name}"
        );
        assert_eq!(case.gencost.as_ref().map(Vec::len), Some(sizes.1), "{name}");
    }
}

#[test]
fn every_matpower_data_file_parses_or_errors_cleanly() {
    // MATPOWER ships 84 data files. Pure data cases must parse; files that compute their data
    //   with code, or aren't cases, must be rejected with an error, never a panic or wrong data
    let Some(case14) = common::matpower_case("case14") else {
        return;
    };
    let data_dir = case14.parent().unwrap();
    let mut parsed = 0;
    for entry in std::fs::read_dir(data_dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().is_some_and(|ext| ext == "m") && Case::from_path(&path).is_ok() {
            parsed += 1;
        }
    }
    assert!(parsed >= 50, "only {parsed} MATPOWER cases parsed");

    // A distribution case that rescales its data with code after the literals
    let rescaled = Case::from_path(data_dir.join("case10ba.m"));
    assert!(parse_message(rescaled)
        .1
        .contains("compute its data with code"));
    // A case with arithmetic in its literals
    let expression = Case::from_path(data_dir.join("case533mt_hi.m"));
    assert!(parse_message(expression).1.contains("`50/3`"));
}

#[test]
fn tiny_case_parses() {
    let case = Tiny::default().parse().unwrap();
    assert_eq!(
        (case.bus.len(), case.generator.len(), case.branch.len()),
        (2, 1, 1)
    );
    assert_eq!(case.gencost, None);
    assert_eq!(case.busnames, None);
    assert_eq!(case.bus_indices()[&2], 1);
    // FromStr is the same as parse
    assert_eq!(Tiny::default().text().parse::<Case>().unwrap(), case);
}

#[test]
fn bus_indices_map_labels_to_positions() {
    let case = Tiny {
        bus: "\t10\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
              \t7\t1\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        gen: "\t10\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
        branch: "\t10\t7\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
        ..Tiny::default()
    }
    .parse()
    .unwrap();
    let indices = case.bus_indices();
    assert_eq!((indices[&10], indices[&7]), (0, 1));
}

#[test]
fn syntax_variations_are_accepted() {
    let text = "%% A comment before the function line\n\
                function mpc = variations % trailing comment\n\
                mpc.version = '2';\n\
                mpc.baseMVA = 100;\n\
                mpc.bus = [ %% comment after the bracket\n\
                \t1, 3, 0, 0, 0, 0, 1, 1, 0, 135, 1, 1.05, 0.95, 99, 99, 99, 99; % extra columns\n\
                \t2 1 50 10 0 0 1 1 0 ...\n\
                \t\t135 1 Inf 0.95 99 99 99 99\n\
                ];\n\
                mpc.gen = [ 1 0 0 100 -100 1 100 1 200 0 0 0 0 0 0 0 0 0 0 0 0 ];\n\
                mpc.branch = [\n\
                \t1\t2\t1e-2\t0.1\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n\
                ];\n\
                mpc.bus_name = { 'Bus 1 % not a comment'; 'Bus ''2''' };\n\
                mpc.gentype = { 'ST' };\n\
                mpc.dcline = [ 1 2 3 ];\n\
                mpc.user_scale = 2.5;\n\
                end\n";
    let case = Case::parse(text).unwrap();
    // Commas, extra columns and a `...` continuation all give two 13-column buses
    assert_eq!(case.bus.len(), 2);
    assert_eq!(case.bus[1].base_kv, 135.0);
    assert_eq!(case.bus[1].vmax_pu, f64::INFINITY);
    assert_eq!(case.branch[0].r, 0.01);
    assert_eq!(
        case.busnames.unwrap(),
        vec!["Bus 1 % not a comment".to_string(), "Bus '2'".to_string()]
    );
}

#[test]
fn statuses_follow_matpower_rules() {
    // A generator is in service when its status is above 0, so -1 is out of service
    let case = Tiny {
        gen: "\t1\t0\t0\t100\t-100\t1\t100\t-1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
        branch: "\t1\t2\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t0\t-360\t360;\n",
        ..Tiny::default()
    }
    .parse()
    .unwrap();
    assert!(!case.generator[0].status);
    assert!(!case.branch[0].status);
}

#[test]
fn gencost_models() {
    let case = Tiny {
        extra: "mpc.gencost = [\n\
                \t1\t0\t0\t3\t0\t0\t50\t1000\t100\t2500;\n\
                \t2\t10\t5\t2\t30\t1\t0\t0\t0\t0;\n\
                ];\n",
        gen: "\t1\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n\
              \t2\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
        ..Tiny::default()
    }
    .parse()
    .unwrap();
    let costs = case.gencost.unwrap();
    assert_eq!(
        costs[0],
        GenCost::PiecewiseLinear {
            startup: 0.0,
            shutdown: 0.0,
            points: vec![(0.0, 0.0), (50.0, 1000.0), (100.0, 2500.0)],
        }
    );
    // A shorter polynomial row, padded with zeros to keep the matrix rectangular
    assert_eq!(
        costs[1],
        GenCost::Polynomial {
            startup: 10.0,
            shutdown: 5.0,
            coeffs: vec![30.0, 1.0],
        }
    );
}

#[test]
fn missing_and_unsupported_fields() {
    let without_version = Tiny::default().text().replace("mpc.version = '2';\n", "");
    assert_eq!(
        Case::parse(&without_version),
        Err(PowerflowError::MissingField { field: "version" })
    );
    let version_1 = Tiny::default().text().replace("'2'", "'1'");
    assert_eq!(
        Case::parse(&version_1),
        Err(PowerflowError::UnsupportedVersion {
            found: "1".to_string()
        })
    );
    let without_bus = Tiny::default()
        .text()
        .replace("mpc.bus = [", "mpc.notbus = [");
    assert_eq!(
        Case::parse(&without_bus),
        Err(PowerflowError::MissingField { field: "bus" })
    );
}

#[test]
fn code_and_expressions_are_rejected() {
    let indexed = Tiny {
        extra: "mpc.bus(:, 3) = mpc.bus(:, 3) / 1e3;\n",
        ..Tiny::default()
    };
    // Lines 1-13 are the function line, version, baseMVA and the three tables
    let (line, message) = parse_message(indexed.parse());
    assert_eq!(line, 14);
    assert!(message.contains("compute its data with code"), "{message}");

    let helper_variable = Tiny {
        extra: "Vbase = 12.66e3;\n",
        ..Tiny::default()
    };
    assert!(parse_message(helper_variable.parse())
        .1
        .contains("unsupported statement"));

    let expression = Tiny::default()
        .text()
        .replace("mpc.baseMVA = 100;", "mpc.baseMVA = 50/3;");
    let (line, message) = parse_message(Case::parse(&expression));
    assert_eq!(line, 3);
    assert!(message.contains("`50/3`"), "{message}");
}

#[test]
fn malformed_tables_report_their_line() {
    // The second bus row is on line 6: function, version, baseMVA, `mpc.bus = [`, first row
    let short_row = Tiny {
        bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
              \t2\t1\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05;\n",
        ..Tiny::default()
    };
    let (line, message) = parse_message(short_row.parse());
    assert_eq!(line, 6);
    assert!(message.contains("first row (line 5)"), "{message}");

    let too_few_columns = Tiny {
        bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05;\n",
        ..Tiny::default()
    };
    let (line, message) = parse_message(too_few_columns.parse());
    assert_eq!(line, 5);
    assert!(message.contains("needs at least 13"), "{message}");

    let fractional_label = Tiny {
        bus: "\t1.5\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        ..Tiny::default()
    };
    assert!(parse_message(fractional_label.parse())
        .1
        .contains("bus number must be a whole number"));

    // The branch table's `];` is the file's last, so without it no `]` ever comes
    let text = Tiny::default().text();
    let unclosed = text.trim_end().strip_suffix("];").unwrap();
    let (line, message) = parse_message(Case::parse(unclosed));
    assert_eq!(line, 11);
    assert!(message.contains("never closes"), "{message}");

    // Without the bus table's `];`, the next table's text is read as bus data
    let run_on = Tiny::default()
        .text()
        .replacen("];\nmpc.gen", "\nmpc.gen", 1);
    assert!(parse_message(Case::parse(&run_on))
        .1
        .contains("`mpc.gen` in a matrix is not a plain number"));

    let assigned_twice = Tiny {
        extra: "mpc.baseMVA = 50;\n",
        ..Tiny::default()
    };
    assert!(parse_message(assigned_twice.parse())
        .1
        .contains("assigned twice (first at line 3)"));
}

#[test]
fn unknown_bus_type() {
    let case = Tiny {
        bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
              \t2\t5\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        ..Tiny::default()
    };
    assert_eq!(
        case.parse(),
        Err(PowerflowError::UnknownBusType { line: 6, code: 5.0 })
    );
}

#[test]
fn validation_errors() {
    let unknown = Tiny {
        branch: "\t1\t99\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
        ..Tiny::default()
    };
    assert_eq!(
        unknown.parse(),
        Err(PowerflowError::UnknownBus {
            table: "branch",
            row: 1,
            label: 99
        })
    );

    let duplicate = Tiny {
        bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
              \t1\t1\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        branch: "\t1\t1\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
        ..Tiny::default()
    };
    assert_eq!(
        duplicate.parse(),
        Err(PowerflowError::DuplicateBus {
            label: 1,
            first_row: 1,
            second_row: 2
        })
    );

    let names = Tiny {
        extra: "mpc.bus_name = { 'only one' };\n",
        ..Tiny::default()
    };
    assert_eq!(
        names.parse(),
        Err(PowerflowError::CountMismatch {
            field: "bus_name",
            expected: "2".to_string(),
            found: 1
        })
    );

    let costs = Tiny {
        extra: "mpc.gencost = [ 2 0 0 2 1 0; 2 0 0 2 1 0; 2 0 0 2 1 0 ];\n",
        ..Tiny::default()
    };
    assert_eq!(
        costs.parse(),
        Err(PowerflowError::CountMismatch {
            field: "gencost",
            expected: "1 or 2".to_string(),
            found: 3
        })
    );

    let zero_base = Tiny::default()
        .text()
        .replace("mpc.baseMVA = 100;", "mpc.baseMVA = 0;");
    assert_eq!(
        Case::parse(&zero_base),
        Err(PowerflowError::InvalidBaseMva { base_mva: 0.0 })
    );
}

#[test]
fn validate_catches_edits() {
    // Edits to a parsed case can be checked again
    let mut case = Tiny::default().parse().unwrap();
    assert_eq!(case.validate(), Ok(()));
    case.branch.push(Branch {
        to_label: 3,
        ..case.branch[0].clone()
    });
    assert_eq!(
        case.validate(),
        Err(PowerflowError::UnknownBus {
            table: "branch",
            row: 2,
            label: 3
        })
    );
}

#[test]
fn errors_display_readably() {
    let err = Case::from_path(Path::new("no/such/case.m")).unwrap_err();
    assert!(err
        .to_string()
        .starts_with("can't read case file no/such/case.m:"));
    let err = PowerflowError::UnknownBus {
        table: "gen",
        row: 3,
        label: 42,
    };
    assert_eq!(
        err.to_string(),
        "gen row 3 refers to bus 42, which doesn't exist"
    );
}

#[test]
fn optional_version_2_columns_default_like_matpower() {
    // A 10-column generator row and an 11-column branch row, as in some version 2 files
    //   (case533mt): MATPOWER's own conversion fills generator columns 11-21 with 0 and the
    //   branch angle limits with -360 and 360
    let case = Tiny {
        gen: "\t1\t0\t0\t100\t-100\t1\t100\t1\t200\t0;\n",
        branch: "\t1\t2\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1;\n",
        ..Tiny::default()
    }
    .parse()
    .unwrap();
    let generator = &case.generator[0];
    assert_eq!((generator.pmax_mw, generator.pmin_mw), (200.0, 0.0));
    assert_eq!((generator.pc1_mw, generator.apf), (0.0, 0.0));
    let branch = &case.branch[0];
    assert_eq!((branch.angmin_deg, branch.angmax_deg), (-360.0, 360.0));

    // Fewer than the columns every version has is still an error
    let short_gen = Tiny {
        gen: "\t1\t0\t0\t100\t-100\t1\t100\t1\t200;\n",
        ..Tiny::default()
    };
    assert!(parse_message(short_gen.parse())
        .1
        .contains("mpc.gen row has 9 columns, but needs at least 10"));
}

#[test]
fn from_bytes_strips_a_utf8_byte_order_mark() {
    let mut bytes = b"\xEF\xBB\xBF".to_vec();
    bytes.extend_from_slice(Tiny::default().text().as_bytes());
    assert_eq!(
        Case::from_bytes(&bytes).unwrap(),
        Tiny::default().parse().unwrap()
    );
}

#[test]
fn from_bytes_falls_back_to_latin1() {
    // A file saved as Latin-1 or Windows-1252, as older MATLAB versions do: "Zürich" with ü
    //   as the single byte 0xFC, which isn't valid UTF-8
    let mut bytes = Tiny {
        extra: "mpc.bus_name = { 'Z",
        ..Tiny::default()
    }
    .text()
    .into_bytes();
    bytes.push(0xFC);
    bytes.extend_from_slice(b"rich'; 'Bern' };\n% caf\xE9 in a comment\n");
    let case = Case::from_bytes(&bytes).unwrap();
    assert_eq!(
        case.busnames.unwrap(),
        vec!["Zürich".to_string(), "Bern".to_string()]
    );

    // Valid UTF-8 is read as UTF-8, not reinterpreted
    let utf8 = Tiny {
        extra: "mpc.bus_name = { 'Zürich'; 'Bern' };\n",
        ..Tiny::default()
    };
    let case = Case::from_bytes(utf8.text().as_bytes()).unwrap();
    assert_eq!(case.busnames.unwrap()[0], "Zürich");
}

#[test]
fn oversized_whole_numbers_are_rejected_not_saturated() {
    // 1e30 is a whole number as an f64, but casting it to usize would silently saturate to
    //   usize::MAX; past 2^53 an f64 can't even represent every whole number
    let huge_label = Tiny {
        bus: "\t1e30\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        ..Tiny::default()
    };
    let (line, message) = parse_message(huge_label.parse());
    assert_eq!(line, 5);
    assert!(
        message.contains("bus number must be a whole number from 1 to"),
        "{message}"
    );

    // 2^53 itself is still accepted
    let largest = Tiny {
        bus: "\t9007199254740992\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n\
              \t2\t1\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95;\n",
        gen: "\t2\t0\t0\t100\t-100\t1\t100\t1\t200\t0;\n",
        branch: "\t9007199254740992\t2\t0.01\t0.1\t0.02\t0\t0\t0\t0\t0\t1;\n",
        ..Tiny::default()
    }
    .parse()
    .unwrap();
    assert_eq!(largest.bus[0].bus_label, 9_007_199_254_740_992);
}

#[test]
fn gencost_n_larger_than_its_row_is_an_error() {
    // An N past the end of the row, including one so large that 4 + 2N would overflow
    for n in ["3", "4503599627370496", "9007199254740992"] {
        let text = format!(
            "{}mpc.gencost = [ 1 0 0 {n} 0 0 ];\n",
            Tiny::default().text()
        );
        let (_, message) = parse_message(Case::parse(&text));
        assert!(
            message.contains(&format!("gencost N is {n}, but the row has only 2 values")),
            "{message}"
        );
    }
}

#[test]
fn nan_in_a_read_column_is_rejected() {
    // A NaN reactance would otherwise only surface as a failed solve
    let case = Tiny {
        branch: "\t1\t2\t0.01\tNaN\t0.02\t0\t0\t0\t0\t0\t1\t-360\t360;\n",
        ..Tiny::default()
    };
    let (line, message) = parse_message(case.parse());
    assert_eq!((line, message.as_str()), (12, "mpc.branch column 4 is NaN"));

    let cost = format!(
        "{}mpc.gencost = [ 2 0 0 2 NaN 0 ];\n",
        Tiny::default().text()
    );
    assert!(parse_message(Case::parse(&cost))
        .1
        .contains("mpc.gencost column 5 is NaN"));

    // Inf is legal (limits use it), and NaN in a column the model doesn't read is ignored
    let ignored = Tiny {
        bus: "\t1\t3\t0\t0\t0\t0\t1\t1\t0\t135\t1\tInf\t0.95\tNaN;\n\
              \t2\t1\t50\t10\t0\t0\t1\t1\t0\t135\t1\t1.05\t0.95\tNaN;\n",
        ..Tiny::default()
    };
    assert_eq!(ignored.parse().unwrap().bus[0].vmax_pu, f64::INFINITY);
}

#[test]
fn a_case_needs_buses() {
    let empty = "function mpc = empty\nmpc.version = '2';\nmpc.baseMVA = 100;\n\
                 mpc.bus = [];\nmpc.gen = [];\nmpc.branch = [];\n";
    assert_eq!(Case::parse(empty), Err(PowerflowError::NoBuses));
}
