# GO Challenge 1 reader checklist

The roadmap for long-term goal 4 in [`LongTermGoals.md`](LongTermGoals.md): reading the input files
of ARPA-E's Grid Optimization (GO) Competition Challenge 1, the SC-ACOPF challenge, into the case
model, so its public datasets can drive goal 3. Each unchecked item in phases G1 to G5 is sized to
be one commit. Phase G0's decisions were made with the maintainer on 2026-11-24.

## Goals

1. **A faithful reader of one instance:** `case.raw`, `case.rop`, `case.inl` and `case.con`,
   with errors that name the file and the line, in the style of the MATPOWER reader.
2. **A validated `Case`** converted from it, which the existing power flows and OPFs solve
   unchanged, plus the data a security-constrained OPF needs that `Case` doesn't hold:
   contingency ratings and voltage limits, switched-shunt ranges, participation factors, and a
   contingency list that maps onto `OpfOutage`.
3. **Out of scope at first:** the later challenges' formats (Challenge 3's single JSON file, with
   multi-period and unit-commitment data), the full PSS/E format beyond what Challenge 1 uses,
   and writing any of these files.

**Conventions carried over from the earlier roadmaps:**
- one commit per checklist item, on `dev/` branches, with the PR's description as its report
- a survey's bar written down before the survey runs
- dataset files fetched by a script, never committed; small results derived from them committed
  with attribution

## What there is to build on

- **The case model.** `Case` (`ybus_powerflow/src/case.rs`) is a faithful copy of a MATPOWER
  file: one load and one shunt per bus, branches as π-models with a tap and phase shift on the
  from side, buses identified by labels. `Case::validate` checks it, and `Case::from_bytes`
  decodes UTF-8 or Latin-1.
- **The MATPOWER reader** (`case/matpower.rs`) reports `PowerflowError::Parse { line, message }`.
  It names no file, since a MATPOWER case is one file.
- **The admittances.** `branch_stamps` and `pi_stamps` in `ac.rs` build each branch's stamps for
  the AC power flow, `Admittances` (shared with the AC OPF) and the OPF derivatives, from the
  same Yf and Yt. The fast-decoupled B′ and B″ are built separately.
- **Outages.** `OpfOutage::branch(row)` and `OpfOutage::generator(row)` (`opf/batch.rs`), and
  `ac_opf_outages` and `dc_opf_outages` over them, take outages by case row.
- **MATPOWER's PSS/E tools.** MATPOWER 8.1, which `scripts/fetch-matpower.sh` fetches, includes
  `psse2mpc`, and Octave runs it. `psse2mpc` folds loads, fixed shunts, switched shunts (at
  their initial value) and branch-level shunts into the buses.

## The format

From the Challenge 1 problem formulation (Appendix A for the files and Appendix C for how they're
read), which OEDI hosts with the datasets (see "The data" below).

- **Syntax (A.1.2, A.1.3).** Fields are separated by commas and strings are in single quotes.
  Each section ends with a line whose first character is `0`, commonly `0 / END OF … DATA`.
  Empty sections are still present, and the multi-section files end with `Q`. Challenge 1
  leaves no field to its default: unused strings are `''` and unused numbers are `0`. `case.con`
  is different: its tokens are separated by spaces, with no quotes.
- **`case.raw` (A.2, C.2 to C.9).** PSS/E revision 33, with two-winding transformers only. The
  sections used, and the fields of each:
  - **case identification:** SBASE
  - **bus:** I, AREA, VM, VA, NVHI, NVLO, EVHI, EVLO
  - **load:** I, ID, STATUS, PL, QL
  - **fixed shunt:** I, ID, STATUS, GL, BL
  - **generator:** I, ID, PG, QG, QT, QB, STAT, PT, PB
  - **non-transformer branch:** I, J, CKT, R, X, B, RATEA, RATEC, ST
  - **transformer**, in four lines: I, J, CKT, MAG1, MAG2, STAT, R12, X12, WINDV1, ANG1, RATA1,
    RATC1, WINDV2
  - **switched shunt:** I, STAT, BINIT, N1 to N8, B1 to B8

  The sections between the transformer and switched shunt sections, such as area interchange
  and the DC lines, aren't used.
- **The transformer model (C.8, equations 68 and 69).** The series admittance is
  g + jb = 1/(R12 + jX12), the ratio τ = WINDV1/WINDV2 and the shift θ = ANG1, with no
  charging. The magnetizing admittance MAG1 + jMAG2 is to ground at bus I and isn't divided by
  τ², so it adds to Yff alone. RATA1 is the base case rating and RATC1 the contingency one.
- **Lines (C.7).** RATEA is the base case rating and RATEC the contingency one. The line-end
  shunts GI, BI, GJ and BJ aren't used.
- **Voltage limits (C.3).** NVHI and NVLO apply in the base case, and EVHI and EVLO in
  contingencies.
- **Switched shunts (C.9).** The susceptance is continuous. Its initial value is BINIT, and its
  bounds are the sums of the negative and of the positive block totals N·B, over the leading
  blocks up to the first zero.
- **Costs, `case.rop` (A.3, C.10).** A generator's (BUS, GENID) record in the generator dispatch
  units gives a DSPTBL. The active power dispatch table with that TBL gives a CTBL, and the
  piecewise linear cost table with that LTBL gives the points: X in MW, Y in $/h.
- **Participation factors, `case.inl` (A.4, C.11).** Field R of each (I, ID) record. It's the
  governor droop in the standard format but a dimensionless participation factor here.
- **Contingencies, `case.con` (A.5, C.12).** Each record is `CONTINGENCY <label>`, then exactly
  one event line, then `END`. The event is either `OPEN BRANCH FROM BUS <i> TO BUS <j> CIRCUIT
  <ckt>`, for a line or a transformer, or `REMOVE UNIT <id> FROM BUS <i>`. A final `END` ends
  the file.

## The data

The datasets are on the Open Energy Data Initiative: Elbert, S., et al. (2024), *ARPA-E Grid
Optimization (GO) Competition Challenge 1* [Data set], OEDI, Pacific Northwest National
Laboratory, <https://doi.org/10.25984/2437761>, under CC BY 4.0. The competition's own site,
gocompetition.energy.gov, was unreachable on 2026-11-24 because of an expired certificate.

- **Layout.** Each bundle has Real-Time ("R") and Offline ("O") parts with the same networks. A
  network folder holds `case.rop`, usually `case.inl`, an `inputfiles.ini` saying where each
  file is, a `readme.txt` and the scenario folders. A scenario folder holds `case.raw` and
  `case.con`, and sometimes its own `case.inl`.
- **The starting point.** Real-Time `case.raw` files carry a base case operating point (VM, VA,
  PG, QG and BINIT) that "should be feasible but not optimal". The Offline files make no such
  claim.
- **The bundles used here** are the two smallest, `Challenge_1_Original_Dataset_1_Scenarios.zip`
  (14 MB) and `Challenge_1_Original_Dataset_2_Scenarios.zip` (10 MB).

## Phase G0: Decisions

- [x] **Decide** how the files map onto the case model. Decided (2026-11-24): grow `Case` now.
      `Branch` gains a magnetizing admittance, so a transformer outage through
      `OpfOutage::branch(row)` takes it out exactly. Folding it into the from bus's shunt, as
      `psse2mpc` does, was the alternative: `Case` would stay unchanged, but an outage would
      leave the shunt behind. This is a breaking change, released as 0.3.0. Everything else
      `Case` can't hold (contingency limits, switched-shunt ranges, element IDs, participation
      factors and contingencies) stays in the reader's own types.
- [x] **Decide** where the reader lives. Decided (2026-11-24): a module of `ybus_powerflow`,
      `case::go`, beside `case::matpower`, which reuses `Case` and `PowerflowError`. A crate of
      its own was the alternative: a leaner power flow crate, but another one to version and
      publish.
- [x] **Decide** how tests get the data. Decided (2026-11-24): fetch it, as MATPOWER's cases
      are fetched. CC BY 4.0 allows committing it with attribution, but fetching keeps the
      repository small and follows the precedent. Hand-written fixtures and small results
      derived from the data are committed, the latter with attribution.

## Phase G1: Magnetizing admittance in `Branch`

- [x] `Branch` gains `g_mag` and `b_mag`, in p.u. on the system base, to ground at the from bus
      and not divided by the tap. The MATPOWER reader sets both to 0, and their docs say
      they're not part of MATPOWER's format. `Case::validate` doesn't check them for NaN, as it
      checks no other value: the readers reject NaN where it's read, and G3's will too. (The
      `Branch` in `tests/data_model_test.rs` is built by struct update, so it needed no change.)
- [x] `branch_stamps` adds `g_mag + j·b_mag` to Yff, after `pi_stamps`'s division by τ², so the
      AC power flow, the AC OPF and the derivatives all include it, from-end flows too. B″ gains
      `-b_mag` at (f, f). The DC models ignore it, as they ignore line charging, and say so.
      Tests:
      - Ybus entries match hand-computed ones.
      - An AC power flow with a magnetizing branch equals the same case with an equal shunt at
        the from bus instead.
      - That transformer's `OpfOutage::branch` equals the case with the branch removed and no
        shunt.

      The MATPOWER reference results are unchanged, since the fields are 0 there
      (`ybus_powerflow/tests/magnetizing_test.rs`).
- [x] `ybus_ffi` reads and writes the new columns, `YBUS_BRANCH_G_MAG` (101) and
      `YBUS_BRANCH_B_MAG` (102), with `ybus.h` regenerated and a C test. Field codes are
      MATPOWER's column numbers, so fields beyond its format take codes from 101 up, clear of
      any column MATPOWER might add. They're constants, not symbols, so `symbols.txt` is
      unchanged. `CHANGELOG.md` records the breaking change under [Unreleased].

## Phase G2: The data

- [x] `scripts/fetch-go-c1.sh` downloads the two bundles into the gitignored `testdata/go-c1/`,
      one directory per bundle, SHA-256 checked, reusing a verified extraction and honoring a
      `GO_C1_DIR` override, as `scripts/fetch-matpower.sh` does. OEDI publishes no checksums,
      so the pins are from the first download (2026-11-26), which a second, separate download
      matched. The README credits the data with the citation above. CI doesn't run it yet:
      decided with the maintainer (2026-11-26), its cached fetch step comes with G4's tests,
      the first that read the data.
- [ ] **Inspect** the data and record it here, before G3 and G4 rely on it:
      - the format of `inputfiles.ini`
      - which bus types (IDE), BASKV, CW, CZ and CM values occur
      - whether any line-end shunt (GI, BI, GJ, BJ) or MBASE is nonzero
      - how many transformers have a nonzero MAG1 or MAG2
      - whether every generator has a cost chain and an `.inl` record
      - each network's size

## Phase G3: The `.raw` reader

- [ ] `PowerflowError::GoParse { file, line, message }`, where `file` is a new `GoFile` (Raw,
      Rop, Inl or Con), displayed as `case.raw line {line}: …`. `ybus_ffi` maps it to
      `YBUS_ERROR_PARSE`. The UTF-8 and Latin-1 decoding moves out of `Case::from_bytes` into a
      helper both readers share.
- [ ] A tokenizer for the comma format (single quotes, `/` starting a comment outside quotes,
      `0` ending a section, `Q` ending the file), whose records keep their line numbers.
- [ ] `RawCase`, with the used fields of each section in file order and the IDs and circuits
      unquoted and trimmed. Unused fields are checked for count but not kept, and unused
      sections are skipped by counting their terminators. Tests:
      - hand-written fixtures in `ybus_powerflow/tests/data/go/` covering every section:
        several loads and shunts on one bus, elements out of service, a transformer with
        magnetizing admittance and an off-nominal ratio and shift, switched-shunt blocks, and
        the skipped sections
      - malformed input, checking the file and line named
      - a `catch_unwind` fuzz test, as `tests/matpower_fuzz_test.rs` does

## Phase G4: Conversion to `Case`

- [ ] `RawCase::to_case` returns a validated `Case` and `GoKeys`. `GoKeys` maps each bus I to
      its row, each generator's (I, ID) to its row, and each branch's (I, J, CKT) to its row,
      with lines before transformers.
      - **Buses:**
        - In-service loads sum into the demand and in-service fixed shunts into the shunt.
        - In-service switched shunts add BINIT to the susceptance.
        - NVHI and NVLO become the voltage limits, and VM and VA are copied.
        - The bus type comes from IDE, or from the rule G2's inspection settles if IDE isn't
          filled.
      - **Generators:**
        - Copied: PG, QG, QT and QB for the reactive limits, PT and PB for the real limits, and
          STAT.
        - `vg_pu` is the bus's VM, since VS isn't used.
        - The capability curve and ramp columns are 0.
      - **Lines:** R, X and B, with RATEA as `rate_a` and RATEC as `rate_c`.
      - **Transformers:** R12 and X12, no charging, WINDV1/WINDV2 as the tap, ANG1 as the
        shift, MAG1 and MAG2 as `g_mag` and `b_mag`, and RATA1 and RATC1 as the ratings.
      - **Angle limits:** ±360.
- [ ] Switched-shunt ranges by C.9's rule, and contingency voltage limits, on `RawCase`. Tests
      check them on the fixtures.
- [ ] **Write the bar** for the base case checks here, before running them on the data: the
      largest bus mismatch, in p.u., allowed at the Real-Time operating point, and how closely
      the AC power flow must reproduce the file's VM and VA.
- [ ] **The base case checks,** on every Real-Time scenario in the two bundles:
      - the AC bus mismatch at the file's operating point, computed with our admittances, which
        checks the transformer and magnetizing model against Challenge 1's own equations
      - an AC power flow, which must converge and reproduce the file's voltages

      Both are judged by the bar.
- [ ] **A cross-check against MATPOWER.** `scripts/generate-go-reference.sh` and a `.m` beside
      `scripts/generate_reference.m` convert each network's base case with `psse2mpc` under
      Octave. The test compares the bus and branch tables, with the magnetizing admittance,
      which `psse2mpc` folds into the buses, compared on its own. The results are committed to
      `ybus_powerflow/tests/data/go-reference/`, with attribution.

## Phase G5: Costs, participation factors, contingencies and the instance

- [ ] **`.rop`:** the three sections used, and each generator's chain resolved to a
      `GenCost::PiecewiseLinear` (startup and shutdown 0), in generator order, into
      `Case.gencost`. A broken chain names the file and line.
- [ ] **`.inl`:** participation factors in generator order. A generator without a record is an
      error or gets 0, as G2's inspection settles.
- [ ] **`.con`:** `Contingency { label, event }`, where `event` is a `ContingencyEvent::Branch
      { from, to, circuit }` or a `ContingencyEvent::Generator { bus, id }`. Through `GoKeys`,
      `Contingency::to_outage` gives the matching `OpfOutage`.
- [ ] **`GoInstance`:** the `RawCase`, the `Case`, the keys, the participation factors and the
      contingencies, read by `from_paths(raw, rop, inl, con)` or by `from_scenario_dir`, which
      finds the files through `inputfiles.ini`.
- [ ] **The survey,** in `docs/benchmarks/go-c1-reader-<date>.md`. It reads every scenario in
      both bundles and records:
      - the reading times
      - the base case checks against the bar
      - G2's counts
      - a DC and an AC OPF of each base case
      - one contingency of each through `ac_opf_outages`
- [ ] **Close the goal:** a README section on reading GO instances, the CHANGELOG entry, and
      goal 4 marked done in `LongTermGoals.md`.
