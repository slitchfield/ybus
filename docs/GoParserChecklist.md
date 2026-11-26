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

### What the data holds (2026-11-26)

Counted over every scenario of both bundles (240 scenarios, 245,702 contingencies) by throwaway
read-only commands, listed at the end of this section. Where the data departs from the spec, the
departure is in bold.

- **Layout.**
  - **Parts and networks:** each bundle has an Offline and a Real-Time part with the same
    networks. Bundle 1 has Networks 01, 03, 05, 07 and 10, with 10 scenarios each. Bundle 2 has
    Networks 01, 03 and 07 with 20 scenarios each, and Network 05 with 10.
  - **`inputfiles.ini`** is an INI file whose one `[INPUTS]` section maps `ROP`, `RAW`, `CON` and
    `INL` to `.` (the network folder) or `x` (each scenario folder). `RAW` and `CON` are always
    `x`. **Network 10 also has `ROP=x`**, so each of its scenarios has its own `case.rop` and
    `case.inl`. Every other network shares both at network level.
  - **A stray `case.raw3`** sits in bundle 2's Real-Time Network_01R-020/scenario_1, beside
    its `case.raw` and differing from it. `inputfiles.ini` doesn't name it, so a reader that
    follows `inputfiles.ini` never sees it.
- **Sizes,** over every scenario of each network, in both parts and both bundles:

  | Network | Scenarios | Buses | Loads | Fixed shunts | Generators | Lines | Transformers | Switched shunts | Contingencies per scenario |
  |---|---|---|---|---|---|---|---|---|---|
  | 01 | 60 | 500 | 200 | 0 | 90 | 468 | 131 | 17 to 19 | 371 to 409 |
  | 03 | 60 | 793 | 568 | 49 | 210 | 769 | 143 | 50 | 86 to 100 |
  | 05 | 40 | 2,000 | 1,125 | 0 to 3 | 544 | 2,359 | 847 | 153 to 160 | 3,143 to 3,187 |
  | 07 | 60 | 2,312 | 1,529 | 121 | 617 | 2,156 | 857 | 201 | 979 to 1,014 |
  | 10 | 20 | 7,977 | 6,671 | 0 to 4,587 | 560 to 570 | 8,896 to 8,900 | 2,548 to 2,693 | 6 to 5,082 | 629 to 2,402 |

  **A network's scenarios aren't all the same network.** Network 10's differ most: 10
  different sets of these counts in its 20 scenarios, matching its per-scenario `case.rop` and
  `case.inl`. Networks 01 and 05 vary only in their shunts, and 03 and 07 not at all.
  (Corrected 2026-11-26: this table first gave one scenario's sizes per network, mislabeled
  as scenario 1's. G3's reader, read over every file, caught it.)
- **Text.**
  - **Encoding:** every file is ASCII, with no byte order mark and no non-ASCII bytes.
  - **Line endings:** CRLF, except that **60 `case.raw` files end lines 2 and 3 (the case
    identification's comment lines) with LF alone**. There are no bare CRs. (`file` reports CR
    endings in 3 files, a false positive of its sampling.)
  - **File endings:** every `case.raw` and `case.rop` ends with a `Q` line. `case.inl` ends with
    a terminator line: `0` without a final newline (8 files), `0 / END OF DATA` (20), or `0`
    followed by a `Q` line without a final newline (8). `case.con` ends with `END`, with or
    without a final newline.
- **`case.raw`'s structure.**
  - **Case identification:** every file is REV 33, with SBASE 100. The comment lines are
    free text.
  - **Sections:**
    - **Terminators:** every terminator starts `0 /` or `0 /` followed directly by text, in
      several wordings, upper or lower case.
    - **Count:** 180 files terminate all 19 sections. **60 end the last, induction machines,
      with `Q` instead of a terminator**, which the spec allows.
    - **Unused sections:** area data is never empty, and zone and owner data sometimes aren't,
      so the skipped sections can hold records.
  - **No data record has a `/` comment.**
  - **Record lengths are uniform:**
    - 13 fields per bus, 14 per load, 5 per fixed shunt, 28 per generator, 24 per line and 26
      per switched shunt
    - transformers in four lines of 21, 3, 17 and 2 fields
  - **Strings:** single-quoted, and IDs and circuits are padded to two characters inside the
    quotes (`'1 '`), so they must be trimmed. Names are empty or blank. Circuits are digits, or
    `TL`.
- **`case.raw`'s values.**
  - **Bus types:**
    - IDE is always filled, and only 1, 2 or 3 occur.
    - Each file has exactly one IDE 3 bus, which always has an in-service generator.
    - **28,344 IDE 2 buses have no in-service generator.** `bustypes` treats them as PQ.
    - BASKV is always positive.
  - **Lines:** GI, BI, GJ and BJ are always 0. RATEA is never 0. R is 0 on 180 lines, and X is
    never 0.
  - **Transformers:**
    - K is always 0, so every transformer is two-winding.
    - **CW, CZ and CM are 1, 1 and 1, except for 5,700 with CZ 2, all in Network 07.** Each of
      those has SBASE12 equal to SBASE (100 MVA), and NOMV1 and NOMV2 equal to their buses'
      BASKV. So R12 and X12 mean the same on either base, and C.8's reading of them is exact.
    - **MAG1 is always 0, and MAG2 is nonzero on 15,180,** which G1 now models.
    - WINDV2 is always 1. WINDV1 is off nominal on 71,048, and ANG1 nonzero on 46.
  - **Generators:** 35,720 of 88,160 are out of service, and MBASE is 0 on 27,160.
  - **Switched shunts:**
    - The blocks have no gaps: no nonzero block follows a zero one.
    - 8 shunts have no blocks, so their range is [0, 0].
    - **690 BINIT values lie outside the C.9 range,** all in the Real-Time parts of Networks 05
      and 07. Of those, 630 are within 0.04 MVAr of the range, 7 are between 0.04 and 1 MVAr
    outside it, and 53 are 1 to 19.9 MVAr outside.
  - **Keys:** no (I, ID) or (I, J, CKT) repeats within a table. No (I, J, CKT) is shared by a
    line and a transformer, and none also appears as (J, I, CKT).
- **Real-Time against Offline.** Matched by key, **the used fields of 50 of the 120 pairs differ
  only in the operating point** (VM, VA, PG, QG and BINIT).
  - **Switched-shunt blocks:** in the other 70, the blocks are split differently (9 × 6 against
    1 × 54 MVAr, for example), usually with the same total. In Network 03, the C.9 ranges of 270
    shunt records, counted over its scenarios, differ slightly (−10 against −9.99 MVAr).
  - **The reference bus:** Network 03's reference bus is a different generator bus in each part
    (29 pairs).

  So each file is consistent on its own, but the parts aren't interchangeable apart from the
  operating point, as the format document implies.
- **`case.rop`.**
  - **Line 1 is data:** its `0` is the modification code record, and the same line is that
    section's terminator. 16 terminated sections follow, then `Q`.
  - **56,980 cost-point lines begin with `0`** (an X of 0 MW), so points must be read by NPAIRS,
    never tested as terminators. That's the spec's warning about `0` lines that are data.
  - **Coverage:** every generator, in service or not, has exactly one dispatch unit and a full
    chain (dispatch unit to dispatch table to cost table), and no record refers to a missing
    one.
  - **Curves:** 2 to 6 points each, with X strictly increasing and Y never decreasing.
- **`case.inl`.** 7 fields per record. Every generator has exactly one record, and none refers to
  a generator that doesn't exist. R is never negative, and it's 0 on 72 in-service generators.
- **`case.con`.**
  - **Events:** each record has at most one event, as the spec says.
    - 207,902 branch outages: 184,506 lines and 23,396 transformers
    - 37,788 generator outages
  - **References:** every event refers to an existing element in the file's own bus order. No
    label repeats within a file. Labels are free-form tokens of letters, digits and punctuation,
    in patterns such as `LINE-n-n-n`, `L_n~n` and `GEN-n-n`, where n stands for a number.
  - **12 records have no event,** all in Network 10's scenario 1 (6 per part).
  - **2 outages name a transformer that's already out of service:** the transformer from bus
    259 to bus 6855, circuit 1, in one Network 10 scenario per part.

<details>
<summary>The commands behind these counts (throwaway, read-only)</summary>

Run from `testdata/go-c1/` after `scripts/fetch-go-c1.sh`, with each script saved outside the data
directories and run as `python3 -I <script> …`, the data passed only as arguments. `raw_values.py`
holds the rough `case.raw` reader the others reuse: they load its `records` function, everything
before its counting. It's a throwaway, not a reference for the Rust reader.

```sh
# Layout: networks and scenarios per part, and inputfiles.ini
for p in */Original_Dataset_*; do echo "$p"; for n in $p/Network_*; do echo "  $(basename $n) $(ls -d $n/scenario_* | wc -l)"; done; done
find . -name inputfiles.ini -exec md5sum {} + | awk '{print $1}' | sort | uniq -c
unzip -Z1 <bundle>.zip | grep -vE '/scenario_[0-9]+/(case\.raw|case\.con)?$' | grep -v '/$'
# Text: encodings, line endings, BOMs, non-ASCII bytes, and how each file type ends
for t in case.raw case.rop case.inl case.con; do find . -name $t -exec file -b {} + | sort | uniq -c; done
find . -name 'case.*' -print0 | xargs -0 python3 -I eol.py
find . -name 'case.*' -print0 | xargs -0 python3 -I bare_cr.py
find . -type f -name 'case.*' -exec grep -lP '[^\x00-\x7F]' {} + | wc -l
for t in case.raw case.rop case.con; do find . -name $t -exec sh -c 'tail -c 12 "$1" | od -An -c' _ {} \; | sort | uniq -c; done
find . -name case.inl -exec sh -c 'tail -n 2 "$1" | tr -d "\r" | tr "\n" "|"; echo' _ {} \; | sort | uniq -c
# case.raw: sections, values, CZ 2, BINIT, Real-Time against Offline
find . -name case.raw -print0 | xargs -0 python3 -I raw_sections.py
find . -name case.raw -path '*scenario_*' -print0 | xargs -0 python3 -I raw_values.py
python3 -I sizes.py raw_values.py   # the size table, over every scenario
find . -name case.raw -path '*scenario_*' -print0 | xargs -0 python3 -I raw_details.py raw_values.py
python3 -I nomv.py raw_values.py $(find . -name case.raw -path '*Network_07*scenario_1/*')
find . -name case.raw -path '*scenario_*' -print0 | xargs -0 python3 -I binit.py raw_values.py
python3 -I ro_compare.py raw_values.py
python3 -I ro_swsh_ide.py raw_values.py
# case.rop and case.inl coverage, then case.con
find . -name case.rop -exec awk '/^ *0 *($|\/)/' {} + | tr -d '\r' | sort | uniq -c
python3 -I rop_inl.py raw_values.py
python3 -I con.py raw_values.py
python3 -I oos_outage.py raw_values.py
for f in $(find . -name case.con); do awk -v f="$f" 'prev ~ /^CONTINGENCY/ && $0 ~ /^ *END *\r?$/ {n++} {prev=$0} END {if (n) print n, f}' "$f"; done
```

`raw_values.py`:

```python
# Rough, throwaway reader of GO Challenge 1 case.raw files (rev 33, two-winding transformers),
#   for counting what the data holds. Not a reference for the Rust reader.
import sys, csv, collections, io, json
SECTIONS = ['bus','load','fshunt','gen','line','xfmr','area','tdc','vsc','icorr','mtdc','msl',
            'zone','iat','owner','facts','swsh','gne','induction']
def records(path):
    lines = open(path, encoding='ascii').read().splitlines()
    out = {s: [] for s in SECTIONS}; sec = 0; i = 3
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith('Q'): break
        if s.startswith('0') and (len(s) == 1 or s[1] in ' /'):
            sec += 1; i += 1; continue
        n = 4 if SECTIONS[sec] == 'xfmr' else 1
        rec = []
        for line in lines[i:i+n]:
            rec.append([f.strip() for f in next(csv.reader([line], quotechar="'", skipinitialspace=True))])
        out[SECTIONS[sec]].append((i+1, rec)); i += n
    return lines[0], out
C = collections.defaultdict(collections.Counter)
sizes = {}
for path in sys.argv[1:]:
    head, r = records(path)
    net = path.split('/')[-3] if 'scenario_' in path else path
    sizes.setdefault(net.replace('O-','?-').replace('R-','?-')[:10], {k: len(v) for k, v in r.items() if k in ('bus','load','fshunt','gen','line','xfmr','swsh')})
    for k in ('bus','load','fshunt','gen','line','swsh'):
        for _, rec in r[k]: C[k+' nfields'][len(rec[0])] += 1
    for _, rec in r['xfmr']: C['xfmr nfields per line'][tuple(len(l) for l in rec)] += 1
    for _, (b,) in r['bus']:
        C['bus IDE'][b[3]] += 1; C['bus BASKV>0'][float(b[2]) > 0] += 1
        C['bus NVHI/NVLO/EVHI/EVLO'][tuple(b[9:13])] += 1
    for _, (l,) in r['load']: C['load STATUS'][l[2]] += 1; C['load ID len/quoted'][(len(l[1]), l[1] != l[1].strip())] += 1
    for _, (f,) in r['fshunt']: C['fshunt STATUS'][f[2]] += 1
    gbus = collections.Counter(); gkeys = collections.Counter()
    for _, (g,) in r['gen']:
        C['gen STAT'][g[14]] += 1; C['gen MBASE==0'][float(g[8]) == 0] += 1
        gkeys[(g[0], g[1])] += 1; gbus[g[0]] += 1
        C['gen ID'][g[1]] += 1
    C['dup gen keys'][sum(1 for v in gkeys.values() if v > 1)] += 1
    ide = {b[0]: b[3] for _, (b,) in r['bus']}
    for bus, ty in ide.items():
        ing = any(g[0] == bus and g[14] == '1' for _, (g,) in r['gen'])
        C['IDE vs in-service gen'][(ty, ing)] += 1
    bkeys = collections.Counter()
    for _, (l,) in r['line']:
        C['line ST'][l[13]] += 1
        C['line GI/BI/GJ/BJ nonzero'][any(float(x) != 0 for x in l[9:13])] += 1
        C['line RATEA==0'][float(l[6]) == 0] += 1; C['line R,X zero'][(float(l[3]) == 0, float(l[4]) == 0)] += 1
        bkeys[(l[0], l[1], l[2])] += 1; C['CKT'][l[2]] += 1
    for _, (a, b, c, d) in r['xfmr']:
        C['xfmr K'][a[2]] += 1; C['xfmr CW,CZ,CM'][(a[4], a[5], a[6])] += 1; C['xfmr STAT'][a[11]] += 1
        C['xfmr MAG1!=0, MAG2!=0'][(float(a[7]) != 0, float(a[8]) != 0)] += 1
        C['xfmr WINDV1!=1, WINDV2!=1, ANG1!=0'][(float(c[0]) != 1, float(d[0]) != 1, float(c[2]) != 0)] += 1
        C['xfmr R12==0, X12==0'][(float(b[0]) == 0, float(b[1]) == 0)] += 1
        C['xfmr NOMV1, NOMV2 nonzero'][(float(c[1]) != 0, float(d[1]) != 0)] += 1
        C['xfmr COD1'][c[6]] += 1
        bkeys[(a[0], a[1], a[3])] += 1; C['CKT'][a[3]] += 1
    lk = {(l[0], l[1], l[2]) for _, (l,) in r['line']}; tk = {(a[0], a[1], a[3]) for _, (a, *_) in r['xfmr']}
    C['dup branch keys within tables'][sum(1 for v in bkeys.values() if v > 1) - len(lk & tk)] += 1
    C['keys shared by a line and a transformer'][len(lk & tk)] += 1
    rev = {(j, i, c) for (i, j, c) in lk | tk}
    C['keys also present reversed'][len(rev & (lk | tk))] += 1
    for _, (s,) in r['swsh']:
        C['swsh STAT'][s[3]] += 1
        blocks = [(int(float(s[10+2*k])), float(s[11+2*k])) for k in range(8)]
        nz = [n * b != 0 for n, b in blocks]
        lead = next((k for k, v in enumerate(nz) if not v), 8)
        C['swsh leading nonzero blocks'][lead] += 1
        C['swsh nonzero block after a zero one'][any(nz[lead:])] += 1
        tot = [n * b for n, b in blocks[:lead]]
        lo = sum(min(0, t) for t in tot); hi = sum(max(0, t) for t in tot); binit = float(s[9])
        C['swsh BINIT within [lo, hi] (1e-6)'][lo - 1e-6 <= binit <= hi + 1e-6] += 1
        C['swsh signs (neg, pos)'][(any(t < 0 for t in tot), any(t > 0 for t in tot))] += 1
    C['unused sections nonempty'][tuple(k for k in ('area','tdc','vsc','icorr','mtdc','msl','zone','iat','owner','facts','gne','induction') if r[k])] += 1
for k, v in C.items():
    items = v.most_common(12)
    print(f'{k}: {items}{" ..." if len(v) > 12 else ""}')
print('sizes:'); [print(' ', k, v) for k, v in sorted(sizes.items())]
```

`sizes.py`:

```python
import sys, collections, glob
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
KEYS = ('bus', 'load', 'fshunt', 'gen', 'line', 'xfmr', 'swsh')
per = collections.defaultdict(list)
for path in sorted(glob.glob('*/Original_Dataset_*/Network_*/scenario_*/case.raw')):
    _, r = ns['records'](path)
    net = [p for p in path.split('/') if p.startswith('Network_')][0]
    scen = path.split('/')[-2]
    per[net[:10]].append((net, scen, tuple(len(r[k]) for k in KEYS)))
for net, rows in sorted(per.items()):
    distinct = sorted(set(sizes for _, _, sizes in rows))
    print(net, len(rows), 'scenarios,', len(distinct), 'distinct size tuples')
    for k, key in enumerate(KEYS):
        vals = [s[k] for _, _, s in rows]
        print(f'   {key:7} {min(vals):>6} to {max(vals):>6}')
    print('   scenario_1 sizes by part:', sorted({(n, s) for n, sc, s in rows if sc == 'scenario_1'}))
```

`eol.py`:

```python
import sys, re, collections
for path in sys.argv[1:]:
    b = open(path, 'rb').read()
    crlf = b.count(b'\r\n'); lf = b.count(b'\n') - crlf; cr = b.count(b'\r') - crlf
    if lf == 0 and cr == 0: continue
    # where: line numbers (by LF) of bare CRs and LF-only lines, and their text
    lines = b.split(b'\n')
    bare_cr = [(i+1, l[:60]) for i, l in enumerate(lines) if b'\r' in l.rstrip(b'\r')]
    lf_only = [i+1 for i, l in enumerate(lines[:-1]) if not l.endswith(b'\r')]
    print(path, f'crlf={crlf} lf_only={lf} bare_cr={cr}', 'lf_only_lines', lf_only[:5], '...' if len(lf_only)>5 else '', 'bare_cr', bare_cr[:3])
```

`bare_cr.py`:

```python
import sys, re
for path in sys.argv[1:]:
    b = open(path, 'rb').read()
    print(path, 'bare CR:', len(re.findall(rb'\r(?!\n)', b)), 'CRCRLF:', b.count(b'\r\r\n'))
```

`raw_sections.py`:

```python
import sys, collections, re
# For each case.raw: the case-ID line, the terminator lines (first char '0' after the 3 ID lines),
#   how many sections, and data lines with '/' outside single quotes.
terms = collections.Counter(); nsec = collections.Counter(); rev = collections.Counter()
slash_data = collections.Counter(); qline = collections.Counter(); first = collections.Counter()
def outside_quotes_slash(line):
    q = False
    for c in line:
        if c == "'": q = not q
        elif c == '/' and not q: return True
    return False
for path in sys.argv[1:]:
    lines = open(path, encoding='ascii').read().splitlines()
    idf = [f.strip() for f in lines[0].split('/')[0].split(',')]
    rev[(len(idf), idf[2] if len(idf) > 2 else None)] += 1
    first[lines[0].split('/')[1].strip()[:30] if '/' in lines[0] else '(no comment)'] += 1
    k = 0; sec = 0
    for line in lines[3:]:
        s = line.strip()
        if s == 'Q' or s.startswith('Q'):
            qline[s[:20]] += 1; break
        if s.startswith('0') and (len(s) == 1 or s[1] in ' /,'):
            # a '0' alone or '0 /...' ends a section (data records start with a bus number > 0)
            terms[re.sub(r'\s+', ' ', s)[:70]] += 1; sec += 1; continue
        if outside_quotes_slash(line): slash_data[sec] += 1
    nsec[sec] += 1
print('case-ID fields, REV:', dict(rev)); print('case-ID comment starts:', dict(first))
print('sections per file:', dict(nsec)); print('Q lines:', dict(qline))
print('data lines with / outside quotes, by section index:', dict(slash_data))
print('terminator lines:'); [print(f'  {n:4} {t}') for t, n in terms.most_common()]
```

`raw_details.py`:

```python
import sys, collections, importlib.util
spec = importlib.util.spec_from_file_location('rv', sys.argv[1]); 
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
records = ns['records']
C = collections.defaultdict(collections.Counter)
for path in sys.argv[2:]:
    head, r = records(path)
    net = [p for p in path.split('/') if p.startswith('Network_')][0]
    sbase = float(head.split('/')[0].split(',')[1])
    C['SBASE'][sbase] += 1
    for _, (a, b, c, d) in r['xfmr']:
        if a[5] == '2':
            C['CZ=2 by network'][net] += 1
            C['CZ=2: SBASE12 vs SBASE'][(float(b[2]), sbase)] += 1
    worst = 0
    for _, (s,) in r['swsh']:
        blocks = [(int(float(s[10+2*k])), float(s[11+2*k])) for k in range(8)]
        lead = next((k for k, (n, bb) in enumerate(blocks) if n * bb == 0), 8)
        tot = [n * bb for n, bb in blocks[:lead]]
        lo = sum(min(0, t) for t in tot); hi = sum(max(0, t) for t in tot); binit = float(s[9])
        if not (lo - 1e-6 <= binit <= hi + 1e-6):
            out = max(lo - binit, binit - hi)
            C['BINIT outside by network'][net] += 1
            C['BINIT outside: STAT'][s[3]] += 1
            worst = max(worst, out); C['BINIT outside, MVAr beyond the range (rounded up)'][int(out) + 1 if out >= 1 else round(out, 3)] += 1
    for _, (l,) in r['load']: pass
for k, v in C.items(): print(f'{k}: {v.most_common(15)}')
```

`nomv.py`:

```python
import sys, collections
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
C = collections.Counter()
for path in sys.argv[2:]:
    _, r = ns['records'](path); kv = {b[0]: float(b[2]) for _, (b,) in r['bus']}
    for _, (a, b, c, d) in r['xfmr']:
        if a[5] == '2':
            n1, n2 = float(c[1]), float(d[1])
            C[('NOMV1==0' if n1 == 0 else 'NOMV1==BASKV_I' if abs(n1 - kv[a[0]]) < 1e-9 else 'NOMV1!=BASKV_I',
               'NOMV2==0' if n2 == 0 else 'NOMV2==BASKV_J' if abs(n2 - kv[a[1]]) < 1e-9 else 'NOMV2!=BASKV_J')] += 1
print(C)
```

`binit.py`:

```python
import sys, collections
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
b = collections.Counter(); worst = 0
for path in sys.argv[2:]:
    _, r = ns['records'](path)
    for _, (s,) in r['swsh']:
        blocks = [(int(float(s[10+2*k])), float(s[11+2*k])) for k in range(8)]
        lead = next((k for k, (n, bb) in enumerate(blocks) if n * bb == 0), 8)
        tot = [n * bb for n, bb in blocks[:lead]]
        lo, hi, x = sum(min(0, t) for t in tot), sum(max(0, t) for t in tot), float(s[9])
        out = max(lo - x, x - hi)
        if out > 1e-6:
            worst = max(worst, out)
            b['<= 0.04 MVAr' if out <= 0.04 + 1e-9 else '0.04 to 1 MVAr' if out < 1 else '>= 1 MVAr'] += 1
print(dict(b), 'total', sum(b.values()), 'worst', round(worst, 3))
```

`ro_compare.py`:

```python
import sys, collections, glob, os
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
# Each table's key, and the fields the operating point may change (bus VM, VA; gen PG, QG; BINIT)
KEYS = {'bus': lambda r: (r[0][0],), 'load': lambda r: (r[0][0], r[0][1]), 'fshunt': lambda r: (r[0][0], r[0][1]),
        'gen': lambda r: (r[0][0], r[0][1]), 'line': lambda r: (r[0][0], r[0][1], r[0][2]),
        'xfmr': lambda r: (r[0][0], r[0][1], r[0][3]), 'swsh': lambda r: (r[0][0],)}
ALLOWED = {('bus', 0): {7, 8}, ('gen', 0): {2, 3}, ('swsh', 0): {9}}
USED = {('bus', 0): {0, 4, 7, 8, 9, 10, 11, 12}, ('load', 0): {0, 1, 2, 5, 6}, ('fshunt', 0): {0, 1, 2, 3, 4},
        ('gen', 0): {0, 1, 2, 3, 4, 5, 14, 16, 17}, ('line', 0): {0, 1, 2, 3, 4, 5, 6, 8, 13},
        ('xfmr', 0): {0, 1, 3, 7, 8, 11}, ('xfmr', 1): {0, 1}, ('xfmr', 2): {0, 2, 3, 5}, ('xfmr', 3): {0},
        ('swsh', 0): {0, 3, 9} | set(range(10, 26))}
def norm(x):
    try: return float(x)
    except ValueError: return x.strip()
C = collections.Counter(); used_diff = collections.Counter(); unused_diff = collections.Counter()
for rpath in sorted(glob.glob('*/Original_Dataset_Real-Time_*/Network_*/scenario_*/case.raw')):
    opath = rpath.replace('Real-Time', 'Offline').replace('R-0', 'O-0').replace('R-1', 'O-1')
    (_, r), (_, o) = ns['records'](rpath), ns['records'](opath)
    clean = True
    for sec, key in KEYS.items():
        rk = {key(rec): rec for _, rec in r[sec]}; ok = {key(rec): rec for _, rec in o[sec]}
        if rk.keys() != ok.keys(): C[f'{sec}: different keys'] += 1; clean = False; continue
        for k, ra in rk.items():
            for li, (la, lb) in enumerate(zip(ra, ok[k])):
                for fi, (x, y) in enumerate(zip(la, lb)):
                    if norm(x) == norm(y): continue
                    if fi in ALLOWED.get((sec, li), set()): C['operating-point field differs'] += 1
                    elif fi in USED.get((sec, li), set()): used_diff[(sec, li, fi)] += 1; clean = False
                    else: unused_diff[(sec, li, fi)] += 1
    C['pairs'] += 1; C['pairs whose used fields differ only in the operating point'] += clean
print(dict(C)); print('used fields differing:', dict(used_diff)); print('unused fields differing:', dict(unused_diff))
```

`ro_swsh_ide.py`:

```python
import sys, collections, glob
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
def rng(s):
    blocks = [(int(float(s[10+2*k])), float(s[11+2*k])) for k in range(8)]
    lead = next((k for k, (n, b) in enumerate(blocks) if n * b == 0), 8)
    tot = [n * b for n, b in blocks[:lead]]
    return round(sum(min(0, t) for t in tot), 6), round(sum(max(0, t) for t in tot), 6), sorted(round(t, 6) for t in tot), blocks[:lead]
C = collections.Counter(); ide = collections.Counter()
for rpath in sorted(glob.glob('*/Original_Dataset_Real-Time_*/Network_*/scenario_*/case.raw')):
    opath = rpath.replace('Real-Time', 'Offline').replace('R-0', 'O-0').replace('R-1', 'O-1')
    (_, r), (_, o) = ns['records'](rpath), ns['records'](opath)
    net = [p for p in rpath.split('/') if p.startswith('Network_')][0]
    ok = {rec[0][0]: rec[0] for _, rec in o['swsh']}
    for _, rec in r['swsh']:
        a, b = rng(rec[0]), rng(ok[rec[0][0]])
        if a[3] != b[3]:
            C[(net, 'same range' if a[:2] == b[:2] else 'different range', 'same block totals' if a[2] == b[2] else 'different totals')] += 1
    ob = {rec[0][0]: rec[0] for _, rec in o['bus']}
    for _, (bus,) in r['bus']:
        if bus[3] != ob[bus[0]][3]: ide[(net, 'R IDE ' + bus[3], 'O IDE ' + ob[bus[0]][3])] += 1
for k, v in sorted(C.items()): print(v, k)
for k, v in ide.items(): print(v, k)
```

`rop_inl.py`:

```python
import sys, csv, collections, glob, os, configparser
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
def fields(line): return [f.strip() for f in next(csv.reader([line], quotechar="'", skipinitialspace=True))]
def rop(path):
    lines = open(path, encoding='ascii').read().splitlines()
    sec = 0; i = 1; out = collections.defaultdict(list); zero_x = 0  # line 1: the modification code
    while i < len(lines):
        s = lines[i].strip()
        if s.startswith('Q'): break
        if s.startswith('0') and (len(s) == 1 or s[1] in ' /'): sec += 1; i += 1; continue
        f = fields(lines[i])
        if sec == 9:  # piecewise linear cost tables: LTBL, LABEL, NPAIRS, then NPAIRS points
            n = int(f[2]); pts = [fields(l) for l in lines[i+1:i+1+n]]
            zero_x += sum(1 for l in lines[i+1:i+1+n] if l.strip().startswith('0'))
            out[sec].append((f, pts)); i += 1 + n; continue
        out[sec].append(f); i += 1
    return out, sec, zero_x
def files(scen):
    net = os.path.dirname(scen); ini = configparser.ConfigParser(); ini.read(os.path.join(net, 'inputfiles.ini'))
    return {t: os.path.join(scen if ini['INPUTS'][t] == 'x' else net, 'case.' + t.lower()) for t in ('RAW', 'ROP', 'INL', 'CON')}
C = collections.Counter(); X = collections.Counter()
for scen in sorted(glob.glob('*/Original_Dataset_*/Network_*/scenario_*')):
    f = files(scen); _, raw = ns['records'](f['RAW']); out, nsec, zx = rop(f['ROP'])
    C[f'rop sections after the code: {nsec}'] += 1; X['PWL point lines starting with 0'] += zx
    gens = {(g[0], g[1]) for _, (g,) in raw['gen']}; ins = {(g[0], g[1]) for _, (g,) in raw['gen'] if g[14] == '1'}
    disp = {(d[0], d[1]): d[3] for d in out[4]}; apd = {a[0]: a[6] for a in out[5]}; pwl = {p[0][0]: p for p in out[9]}
    C['dispatch units == generators'] += set(disp) == gens
    chain = {g: pwl.get(apd.get(disp.get(g))) for g in gens}
    X['generators without a full chain'] += sum(1 for v in chain.values() if v is None)
    X['in-service generators without a full chain'] += sum(1 for g in ins if chain[g] is None)
    X['dispatch units for unknown generators'] += len(set(disp) - gens)
    for g, p in chain.items():
        if p is None: continue
        xs = [float(a) for a, _ in p[1]]; ys = [float(b) for _, b in p[1]]
        X['cost curves'] += 1; X['NPAIRS'] += 0
        C[f'NPAIRS {len(xs)}'] += 1
        X['X not increasing'] += any(b <= a for a, b in zip(xs, xs[1:]))
        X['Y decreasing'] += any(b < a for a, b in zip(ys, ys[1:]))
    # INL: I, ID, H, PMAX, PMIN, R, D, ended by 0
    inl = {}
    for line in open(f['INL'], encoding='ascii').read().splitlines():
        s = line.strip()
        if not s or s.startswith('0') or s.startswith('Q'): break
        v = fields(line); inl[(v[0], v[1])] = float(v[5]); C[f'inl fields {len(v)}'] += 1
    X['generators without an inl record'] += len(gens - set(inl)); X['in-service without inl'] += len(ins - set(inl))
    X['inl records for unknown generators'] += len(set(inl) - gens)
    X['inl R == 0 (in service)'] += sum(1 for g in ins if inl.get(g) == 0.0); X['inl R < 0'] += sum(1 for v in inl.values() if v < 0)
    C['scenarios'] += 1
print(dict(C)); print(dict(X))
```

`con.py`:

```python
import sys, collections, glob, os, configparser, re
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
def files(scen):
    net = os.path.dirname(scen); ini = configparser.ConfigParser(); ini.read(os.path.join(net, 'inputfiles.ini'))
    return {t: os.path.join(scen if ini['INPUTS'][t] == 'x' else net, 'case.' + t.lower()) for t in ('RAW', 'CON')}
C = collections.Counter(); shapes = collections.Counter(); labels = collections.Counter(); per = {}
for scen in sorted(glob.glob('*/Original_Dataset_*/Network_*/scenario_*')):
    f = files(scen); _, raw = ns['records'](f['RAW'])
    lines = {(l[0], l[1], l[2]): l[13] for _, (l,) in raw['line']}
    xf = {(a[0], a[1], a[3]): a[11] for _, (a, *_) in raw['xfmr']}
    gens = {(g[1], g[0]): g[14] for _, (g,) in raw['gen']}
    toks = [l.split() for l in open(f['CON'], encoding='ascii').read().splitlines()]
    toks = [t for t in toks if t]
    i = 0; n = 0; seen = collections.Counter()
    while i < len(toks):
        t = toks[i]
        if t == ['END']:
            C['final END'] += 1; i += 1; C['lines after final END'] += len(toks) - i; break
        assert t[0] == 'CONTINGENCY', (scen, t)
        label = t[1]; seen[label] += 1; labels[re.sub(r'\d+', 'N', label)] += 1; C[f'start line tokens {len(t)}'] += 1
        body = []; i += 1
        while toks[i] != ['END']: body.append(toks[i]); i += 1
        i += 1; n += 1; C[f'events per record {len(body)}'] += 1
        for e in body:
            shapes[' '.join(w if w.isalpha() and w.isupper() else 'x' for w in e)] += 1
            if e[:2] == ['OPEN', 'BRANCH']:
                k = (e[4], e[7], e[9])
                if k in lines: C['branch event -> line'] += 1; C[f'  line in service {lines[k]}'] += 1
                elif k in xf: C['branch event -> transformer'] += 1; C[f'  xfmr in service {xf[k]}'] += 1
                elif (e[7], e[4], e[9]) in lines or (e[7], e[4], e[9]) in xf: C['branch event reversed'] += 1
                else: C['branch event -> unknown'] += 1
            elif e[:2] == ['REMOVE', 'UNIT']:
                k = (e[2], e[5])
                C['generator event' if k in gens else 'generator event -> unknown'] += 1
                if k in gens: C[f'  gen in service {gens[k]}'] += 1
    C['duplicate labels'] += sum(1 for v in seen.values() if v > 1)
    net = [p for p in scen.split('/') if p.startswith('Network_')][0][:10]; per.setdefault(net, set()).add(n)
print(dict(C)); print('event shapes:', dict(shapes)); print('label shapes:', labels.most_common(8)); print('contingencies per scenario, by network:', {k: (min(v), max(v)) for k, v in sorted(per.items())})
```

`oos_outage.py`:

```python
import sys, glob, os, configparser, collections
src = open(sys.argv[1]).read().split('C = collections.defaultdict')[0]; ns = {}; exec(src, ns)
for scen in sorted(glob.glob('*/Original_Dataset_*/Network_*/scenario_*')):
    net = os.path.dirname(scen); ini = configparser.ConfigParser(); ini.read(os.path.join(net, 'inputfiles.ini'))
    con = os.path.join(scen if ini['INPUTS']['CON'] == 'x' else net, 'case.con')
    _, raw = ns['records'](os.path.join(scen, 'case.raw'))
    off = {(a[0], a[1], a[3]) for _, (a, *_) in raw['xfmr'] if a[11] == '0'} | {(l[0], l[1], l[2]) for _, (l,) in raw['line'] if l[13] == '0'}
    offg = {(g[1], g[0]) for _, (g,) in raw['gen'] if g[14] == '0'}
    for line in open(con, encoding='ascii'):
        t = line.split()
        if t[:2] == ['OPEN', 'BRANCH'] and (t[4], t[7], t[9]) in off: print('out-of-service branch outage:', scen, t[4], t[7], t[9])
        if t[:2] == ['REMOVE', 'UNIT'] and (t[2], t[5]) in offg: print('out-of-service generator outage:', scen)
```

</details>

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
- [x] **Inspect** the data and record it, before G3 and G4 rely on it: "What the data holds",
      under "The data" above. Decided with the maintainer (2026-11-26), from those findings:
      - **Bus types: copy IDE.** It's always filled, with one reference bus per file that has
        an in-service generator. `bustypes` (`topology.rs`) already treats a PV bus with no
        in-service generator as PQ, as MATPOWER does, so no rule of the reader's own is needed.
        Deriving types from the generators was the alternative, but it would ignore data every
        file has.
      - **A generator without a `case.inl` record is an error,** naming the file and the
        generator. It never happens in the data, so a missing record means a broken or
        mismatched file. Reading it as a participation factor of 0 was the alternative.
      - **CZ other than 1 is read only where it's equivalent to 1:** CZ 2 with SBASE12 equal to
        SBASE and NOMV1 and NOMV2 each 0 or their bus's BASKV, as all 5,700 in the data are.
        Anything else is an error naming the line, rather than a misreading. Converting by
        PSS/E's rules would go beyond the formulation, with nothing in the data to test it.
      - **BINIT is kept as given,** outside its C.9 range or not, and the range is recorded
        beside it. The operating point is a starting point, and G4's mismatch check compares
        against the file's. Clamping it into the range was the alternative.

## Phase G3: The `.raw` reader

- [x] `PowerflowError::GoParse { file, line, message }`, where `file` is a new `GoFile` (Raw,
      Rop, Inl or Con), displayed as `case.raw line {line}: …`. The UTF-8 and Latin-1 decoding
      moves out of `Case::from_bytes` into `decode`, a helper both readers share. `ybus_ffi` is
      unchanged, as decided with the maintainer (2026-11-26). The C library can't read GO files
      yet, so `GoParse` falls to `YBUS_ERROR_OTHER`, as `status.rs`'s test expects of a variant
      added later, until the C library gains GO entry points. Mapping it to `YBUS_ERROR_PARSE`,
      as first planned, would break that test's rule of one code per variant.
- [x] A tokenizer for the comma format (single quotes, `/` starting a comment outside quotes,
      `0` ending a section, `Q` ending the file), whose records keep their line numbers
      (`case/go/text.rs`, with unit tests). Strict, as decided with the maintainer (2026-11-26):
      what A.1.3 rules out, such as fields separated by blanks, double quotes and empty fields
      left to defaults, is an error. From G2's findings:
      - Lines end with CRLF or LF, mixed within a file.
      - The case identification's lines 2 and 3 are free text.
      - A terminator is `0` alone or followed by a `/` comment, with or without a space.
      - `Q` may stand in for the terminators still owed, as when the induction machine section
        isn't terminated.
- [x] `RawCase`, with the used fields of each section in file order and the IDs and circuits
      unquoted and trimmed (the data pads them, as in `'1 '`). Unused fields are checked for
      count but not kept, and unused sections are skipped by counting their terminators.
      They can hold records: area data always does.
      - **The types:** `RawCase` and its record types are `#[non_exhaustive]`, as decided with
        the maintainer (2026-11-26), and each record keeps the line it came from, for G4's
        errors (`case/go/raw.rs`).
      - **Two checks moved here from G4,** since they're about single fields: IDE must be 1 to
        4, and K must be 0, because a three-winding record is five lines.
      - **The rest stays for G4:** the CW, CZ and CM rules, and what a nonzero line-end shunt
        means (the data has none). So do all cross-record checks, such as duplicate keys and
        unknown buses.

      Tests:
      - hand-written fixtures in `ybus_powerflow/tests/data/go/` covering every section:
        several loads and shunts on one bus, elements out of service, a transformer with
        magnetizing admittance and an off-nominal ratio and shift, switched-shunt blocks, and
        skipped sections with records in them. They also cover the text G2 found: mixed line
        endings, padded IDs, both terminator spellings, and `Q` in place of the last
        terminators. The fixture is committed with LF endings, and `tests/go_raw_test.rs` builds
        the CRLF and mixed variants in code, so git's line-ending handling can't change it.
      - malformed input, checking the file and line named (`tests/go_raw_syntax_test.rs`, 23
        edits of the fixture)
      - a `catch_unwind` fuzz test, as `tests/matpower_fuzz_test.rs` does
        (`tests/go_raw_fuzz_test.rs`)

## Phase G4: Conversion to `Case`

- [ ] `RawCase::to_case` returns a validated `Case` and `GoKeys`. `GoKeys` maps each bus I to
      its row, each generator's (I, ID) to its row, and each branch's (I, J, CKT) to its row,
      with lines before transformers.
      - **Buses:**
        - In-service loads sum into the demand and in-service fixed shunts into the shunt.
        - In-service switched shunts add BINIT to the susceptance, as given, even outside its
          C.9 range (decided in G2).
        - NVHI and NVLO become the voltage limits, and VM and VA are copied.
        - The bus type is IDE's, copied (decided in G2). G3's reader already refuses an IDE
          other than 1 to 4.
      - **Generators:**
        - Copied: PG, QG, QT and QB for the reactive limits, PT and PB for the real limits, and
          STAT.
        - `vg_pu` is the bus's VM, since VS isn't used.
        - The capability curve and ramp columns are 0.
      - **Lines:** R, X and B, with RATEA as `rate_a` and RATEC as `rate_c`.
      - **Transformers:** R12 and X12, no charging, WINDV1/WINDV2 as the tap, ANG1 as the
        shift, MAG1 and MAG2 as `g_mag` and `b_mag`, and RATA1 and RATC1 as the ratings.
        CW and CM must be 1, and CZ 1, or 2 where that's equivalent (SBASE12 equal to SBASE,
        and NOMV1 and NOMV2 each 0 or their bus's BASKV, decided in G2). Anything else is an
        error naming the line. G3's reader already refuses a nonzero K.
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
      `Case.gencost`. A broken chain names the file and line. From G2's findings:
      - Line 1 is the modification code record, and is also its section's terminator.
      - Cost points are read by NPAIRS, never tested as terminators, since 56,980 begin with
        `0`.
- [ ] **`.inl`:** participation factors in generator order. A generator without a record is an
      error naming the file and the generator (decided in G2). The file ends with a terminator,
      or with `0` and then `Q`, with or without a final newline.
- [ ] **`.con`:** `Contingency { label, event }`, where `event` is a `ContingencyEvent::Branch
      { from, to, circuit }` or a `ContingencyEvent::Generator { bus, id }`. Through `GoKeys`,
      `Contingency::to_outage` gives the matching `OpfOutage`. From G2's findings:
      - **A record can have no event** (12 in the data), so `event` is an `Option`, and
        `to_outage` of one gives `OpfOutage::default()`, an outage of nothing: the base case.
      - **An outage can name an element already out of service** (2 in the data). It's kept as
        written, and `to_outage` maps it to the element's row like any other. `OpfOutage`
        already defines that case: "a row already out of service stays out".
      - **A label** is any token without spaces.
- [ ] **`GoInstance`:** the `RawCase`, the `Case`, the keys, the participation factors and the
      contingencies, read by `from_paths(raw, rop, inl, con)` or by `from_scenario_dir`. That
      finds the files through the network folder's `inputfiles.ini`, whose `[INPUTS]` section
      maps each file type to `.` (the network folder) or `x` (the scenario folder). It ignores
      files that `inputfiles.ini` doesn't name, such as the stray `case.raw3`.
- [ ] **The survey,** in `docs/benchmarks/go-c1-reader-<date>.md`. It reads every scenario in
      both bundles and records:
      - the reading times
      - the base case checks against the bar
      - G2's counts
      - a DC and an AC OPF of each base case
      - one contingency of each through `ac_opf_outages`
- [ ] **Close the goal:** a README section on reading GO instances, the CHANGELOG entry, and
      goal 4 marked done in `LongTermGoals.md`.
