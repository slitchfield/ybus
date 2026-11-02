% Generates DC and AC power flow reference results by running MATPOWER's rundcpf and runpf
%   (with its default options: Newton's method, tolerance 1e-8, reactive limits not enforced;
%   then again with fast-decoupled XB and BX). Run
%   through
%   scripts/generate-reference.sh, which fetches MATPOWER and passes its directory, the output
%   directory and any cases as arguments.
%
% By default, generates the committed results in powerflow/tests/data/matpower-reference/,
%   which are only for cases derived from public sources (see docs/PowerFlowChecklist.md,
%   Phase 1.5). Cases named as arguments are for local comparisons, and their results are never
%   committed. A case MATPOWER fails on is reported and skipped, and the script ends with an
%   error after writing the others.
%
% For each case, writes:
%   <case>_dcpf_bus.csv     bus number and voltage angle (degrees)
%   <case>_dcpf_branch.csv  branch row (1-based, in case-file order), from and to bus, status,
%                           and real power flow at the from end (MW)
%   <case>_dcpf_gen.csv     generator row (1-based), bus, and real power output (MW)
%   <case>_acpf_bus.csv     bus number, voltage magnitude (p.u.) and angle (degrees)
%   <case>_acpf_branch.csv  branch row, from and to bus, status, and real and reactive power
%                           flows at the from end and at the to end (MW, MVAr)
%   <case>_acpf_gen.csv     generator row, bus, and real and reactive power output (MW, MVAr)
%   <case>_acpf_summary.csv Newton iterations and whether runpf converged
%   <case>_qlim1_*.csv, <case>_qlim2_*.csv
%                           the same for Newton's method with reactive limits enforced
%                           (pf.enforce_q_lims 1 and 2, through the legacy runpf), plus
%                           <case>_qlimN_limited.csv, the generator rows held at a limit
%   <case>_fdxb_*.csv, <case>_fdbx_*.csv
%                           the same for runpf's fast-decoupled XB and BX methods (pf.alg
%                           FDXB and FDBX, through the legacy fdpf rather than MP-Core), whose
%                           iterations are P half-iterations
% Values are printed with 17 significant digits, enough to round-trip a double exactly.

% Marks this file as a script, so Octave accepts the function definition that follows. Octave
%   needs a script's functions defined before they are called.
1;

function write_csv(path, header, row_format, data)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', header);
    % With no rows, fprintf would still print the format once, leaving a blank line
    if ~isempty(data)
        fprintf(fid, [row_format '\n'], data');
    end
    fclose(fid);
end

args = argv();
if numel(args) < 2
    error('usage: generate_reference.m <matpower-dir> <output-dir> [<case>...]');
end
[matpower_dir, out_dir] = deal(args{1}, args{2});
addpath(genpath(matpower_dir));
define_constants;

if numel(args) > 2
    cases = args(3:end);
else
    cases = {'case14', 'case30', 'case118'};
end
mpopt = mpoption('verbose', 0, 'out.all', 0);

failed = {};
for k = 1:numel(cases)
    name = cases{k};
    try
        [result, success] = rundcpf(name, mpopt);
    catch err
        fprintf(2, '%s: %s\n', name, err.message);
        success = false;
    end
    if ~success
        fprintf(2, 'rundcpf failed on %s\n', name);
        failed{end + 1} = name;
        continue;
    end

    write_csv(fullfile(out_dir, [name '_dcpf_bus.csv']), 'bus,va_deg', '%d,%.17g', ...
        [result.bus(:, BUS_I), result.bus(:, VA)]);
    rows = (1:size(result.branch, 1))';
    write_csv(fullfile(out_dir, [name '_dcpf_branch.csv']), 'row,from,to,status,pf_mw', ...
        '%d,%d,%d,%d,%.17g', [rows, result.branch(:, [F_BUS, T_BUS, BR_STATUS, PF])]);
    rows = (1:size(result.gen, 1))';
    write_csv(fullfile(out_dir, [name '_dcpf_gen.csv']), 'row,bus,pg_mw', '%d,%d,%.17g', ...
        [rows, result.gen(:, [GEN_BUS, PG])]);
    printf('%s: %d buses, %d branches, %d generators\n', name, size(result.bus, 1), ...
        size(result.branch, 1), size(result.gen, 1));

    % The AC power flow: runpf's defaults (Newton's method), then fast-decoupled XB and BX
    %   Fast-decoupled runs bypass MP-Core (exp.use_legacy_core), which MATPOWER 8.1 otherwise
    %   uses under Octave too: its fast-decoupled results differ from the legacy fdpf's by up
    %   to 2.3e-7 (case14, FDXB) at the same iteration counts, and powerflow follows fdpf. For
    %   Newton's method the two agree within 6.4e-12, so it keeps runpf's defaults.
    %   Runs with reactive limits enforced (pf.enforce_q_lims 1 and 2) also use the legacy core:
    %   over all of MATPOWER's cases, MP-Core's results differ from it (by up to 71 MW on
    %   case_SyntheticUSA, and in convergence on case6470rte), though on case14, case30 and
    %   case118 the two agree within 5.1e-12.
    for run = {{'acpf', 'NR', 0, 0}, {'fdxb', 'FDXB', 1, 0}, {'fdbx', 'FDBX', 1, 0}, ...
               {'qlim1', 'NR', 1, 1}, {'qlim2', 'NR', 1, 2}}
        [prefix, alg, legacy, qlim] = deal(run{1}{1}, run{1}{2}, run{1}{3}, run{1}{4});
        try
            [result, success] = runpf(name, mpoption(mpopt, 'pf.alg', alg, ...
                'exp.use_legacy_core', legacy, 'pf.enforce_q_lims', qlim));
        catch err
            fprintf(2, '%s (%s): %s\n', name, alg, err.message);
            success = false;
        end
        if ~success
            fprintf(2, 'runpf (%s) failed on %s\n', alg, name);
            failed{end + 1} = [name ' (' alg ')'];
            continue;
        end
        write_csv(fullfile(out_dir, [name '_' prefix '_bus.csv']), 'bus,vm_pu,va_deg', ...
            '%d,%.17g,%.17g', [result.bus(:, BUS_I), result.bus(:, VM), result.bus(:, VA)]);
        rows = (1:size(result.branch, 1))';
        write_csv(fullfile(out_dir, [name '_' prefix '_branch.csv']), ...
            'row,from,to,status,pf_mw,qf_mvar,pt_mw,qt_mvar', ...
            '%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g', ...
            [rows, result.branch(:, [F_BUS, T_BUS, BR_STATUS, PF, QF, PT, QT])]);
        rows = (1:size(result.gen, 1))';
        write_csv(fullfile(out_dir, [name '_' prefix '_gen.csv']), 'row,bus,pg_mw,qg_mvar', ...
            '%d,%d,%.17g,%.17g', [rows, result.gen(:, [GEN_BUS, PG, QG])]);
        if qlim
            % The generators held at a limit, which runpf doesn't report: those in service at a bus
            %   runpf made PQ (it was REF or PV in the case) whose output is exactly at one of its
            %   limits. A generator whose limits are equal is at them whenever its bus is PQ, so
            %   it counts only if no generator with unequal limits at its bus does: otherwise that
            %   one made the bus PQ (in case2736sp, case2746wp and case2746wop)
            mpc = loadcase(name);
            [~, bus_row] = ismember(result.gen(:, GEN_BUS), result.bus(:, BUS_I));
            was = mpc.bus(bus_row, BUS_TYPE);
            now = result.bus(bus_row, BUS_TYPE);
            at_limit = result.gen(:, GEN_STATUS) > 0 & (was == REF | was == PV) & now == PQ ...
                & (result.gen(:, QG) == result.gen(:, QMAX) | result.gen(:, QG) == result.gen(:, QMIN));
            fixed = result.gen(:, QMAX) == result.gen(:, QMIN);
            varied_buses = result.gen(at_limit & ~fixed, GEN_BUS);
            held = find(at_limit & ~(fixed & ismember(result.gen(:, GEN_BUS), varied_buses)));
            write_csv(fullfile(out_dir, [name '_' prefix '_limited.csv']), 'row', '%d', held);
        end
        % The summary last, so that its presence marks a complete run
        write_csv(fullfile(out_dir, [name '_' prefix '_summary.csv']), 'iterations,success', ...
            '%d,%d', [result.iterations, success]);
        printf('%s: %s (enforce_q_lims %d) converged in %d iterations\n', name, alg, qlim, ...
            result.iterations);
    end
end

% The run's provenance, so a regeneration with a different MATPOWER is visible in the diff
fid = fopen(fullfile(out_dir, 'VERSIONS'), 'w');
fprintf(fid, 'MATPOWER %s\nOctave %s\n', mpver, OCTAVE_VERSION);
fclose(fid);

if ~isempty(failed)
    error('MATPOWER failed on %d runs over %d cases: %s', numel(failed), numel(cases), ...
        strjoin(failed, ', '));
end
