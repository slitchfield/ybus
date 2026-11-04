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
%   <case>_islands.csv      bus number and its island from find_islands (1-based, in
%                           MATPOWER's order, largest first), or 0 for a bus with no in-service
%                           branch
%   <case>_ptdf.csv         makePTDF's PTDF with the default slack (the first REF bus), after
%                           ext2int: a header of the columns' bus numbers, then one line per
%                           in-service branch, its 1-based row first. Every bus on cases of up
%                           to 300 buses, otherwise an evenly spaced sample of 50 (bus_idx)
%   <case>_ptdf_distributed.csv
%                           the same with the slack distributed in proportion to each bus's
%                           Pd (makePTDF's vector slack), where some bus has load; for a sample,
%                           each column is the transfer from its bus to the weights (txfr)
%   <case>_lodf.csv         makeLODF's LODF: a header of the outaged branches' rows, then one
%                           line per in-service branch, its row first. Every branch up to 300
%                           buses, otherwise 50 evenly spaced, from makePTDF's transfer columns
%                           (txfr) with makeLODF's formula. The committed results have these
%                           three for case14 and case30 only, and none for a case with
%                           several islands, which a single slack can't balance
%   <case>_bridges.csv      branch rows (1-based, increasing) find_bridges reports as bridges,
%                           only where it succeeds: it fails on cases with parallel branches
%                           whose summed rows pass the end of the table, with more buses than
%                           branches, or deep enough to pass Octave's recursion limit
%   <case>_dcopf_bus.csv    the DC optimal power flow (rundcopf, with opf.dc.solver MIPS forced,
%                           so that MIPS's failures are recorded too): bus number, voltage angle
%                           (degrees) and real power price LAM_P ($/MWh)
%   <case>_dcopf_gen.csv    generator row, bus, status, real power output (MW), and the
%                           multipliers of its lower and upper real power limits ($/MWh)
%   <case>_dcopf_branch.csv branch row, from and to bus, status, real power flow at the from end
%                           (MW), the multipliers of its flow limit at the from and to ends
%                           ($/MWh), and of its angle difference limits ($/degree)
%   <case>_dcopf_summary.csv
%                           the objective ($/hr), MIPS's iterations, whether it succeeded, and
%                           the objective and solver of rundcopf with its default solver (GLPK
%                           for linear costs), as a check of the optimum
%   <case>_dcopf_history.csv
%                           MIPS's termination measures, barrier parameter, step size, objective
%                           and step lengths at the start and after each iteration
%   <case>_opf_*.csv        the same for the AC optimal power flow (runopf, MIPS by default),
%                           with voltage magnitudes (p.u.), reactive power and its price LAM_Q,
%                           the multipliers of the voltage magnitude and reactive power limits,
%                           and flows at both ends; the summary has no default-solver columns
% The OPF runs use MATPOWER's defaults (MP-Core), which agree with the legacy code on OPF.
%   Their results are written even when MIPS fails, since parity includes its failures.
% Values are printed with 17 significant digits, enough to round-trip a double exactly.
%
% YBUS_REFERENCE_RUNS, if set, chooses which of the runs islands, ptdf, dcpf, acpf, fdxb,
%   fdbx, qlim1, qlim2, dcopf and opf to make, separated by commas (default: all), so that new results don't
%   need the slow ones again. The AC runs follow dcpf only if it's chosen and succeeds.

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

% A dense matrix with a header of column labels, and each line led by its row label
function write_matrix(path, column_labels, row_labels, data)
    fid = fopen(path, 'w');
    fprintf(fid, 'row');
    fprintf(fid, ',%d', column_labels);
    fprintf(fid, '\n');
    format = ['%d' repmat(',%.17g', 1, size(data, 2)) '\n'];
    fprintf(fid, format, [row_labels(:), full(data)]');
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
    % Local results: PTDF and LODF for every case named
    ptdf_cases = cases;
else
    cases = {'case14', 'case30', 'case118'};
    % Committed results: case118's PTDF and LODF would take 1.6 MB, about ten times the
    %   packed repository, so they're compared locally instead
    ptdf_cases = {'case14', 'case30'};
end
mpopt = mpoption('verbose', 0, 'out.all', 0);
chosen = getenv('YBUS_REFERENCE_RUNS');
if isempty(chosen)
    chosen = {'islands', 'ptdf', 'dcpf', 'acpf', 'fdxb', 'fdbx', 'qlim1', 'qlim2', 'dcopf', ...
        'opf'};
else
    chosen = strsplit(chosen, ',');
end
wanted = @(run) any(strcmp(run, chosen));

failed = {};
for k = 1:numel(cases)
    name = cases{k};

    % Islands and bridges, which need no solve
    if wanted('islands')
        mpc = loadcase(name);
        [groups, ~] = find_islands(mpc);
        island = zeros(size(mpc.bus, 1), 1);
        for g = 1:numel(groups)
            island(groups{g}) = g;
        end
        write_csv(fullfile(out_dir, [name '_islands.csv']), 'bus,island', '%d,%d', ...
            [mpc.bus(:, BUS_I), island]);
        try
            [~, bridges] = find_bridges(mpc);
            write_csv(fullfile(out_dir, [name '_bridges.csv']), 'row', '%d', ...
                sort(vertcat(bridges{:})));
        catch err
            fprintf(2, '%s: find_bridges failed: %s\n', name, err.message);
        end
    end

    % PTDF and LODF, in ext2int's internal numbering, labelled with the case's bus numbers
    %   and 1-based branch rows. The full matrices are dense, so larger cases get a sample of
    %   columns, computed column by column as makePTDF's bus_idx and txfr do
    if wanted('ptdf') && any(strcmp(name, ptdf_cases))
        mpc = ext2int(loadcase(name));
    end
    % A single slack can't balance several islands: makePTDF's matrix is then singular
    if wanted('ptdf') && any(strcmp(name, ptdf_cases)) && numel(find_islands(mpc)) > 1
        fprintf(2, '%s: %d islands, so no PTDF or LODF\n', name, numel(find_islands(mpc)));
        for file = {'_ptdf.csv', '_ptdf_distributed.csv', '_lodf.csv'}
            if exist(fullfile(out_dir, [name file{1}]), 'file')
                delete(fullfile(out_dir, [name file{1}]));
            end
        end
    elseif wanted('ptdf') && any(strcmp(name, ptdf_cases))
        nb = size(mpc.bus, 1);
        nl = size(mpc.branch, 1);
        labels = mpc.order.bus.i2e;
        rows = mpc.order.branch.status.on;
        slack = find(mpc.bus(:, BUS_TYPE) == REF, 1);
        if nb <= 300
            buses = (1:nb)';
            outages = (1:nl)';
        else
            buses = unique(round(linspace(1, nb, 50)))';
            outages = unique(round(linspace(1, nl, 50)))';
        end
        if nb <= 300
            H = makePTDF(mpc, slack);
        else
            H = makePTDF(mpc, slack, buses);
        end
        write_matrix(fullfile(out_dir, [name '_ptdf.csv']), labels(buses), rows, H);
        w = mpc.bus(:, PD);
        if any(w)
            if nb <= 300
                Hd = makePTDF(mpc, w);
            else
                % makePTDF's vector slack with bus_idx adds a column for every bus with
                %   weight, which on the largest cases is tens of GB. A distributed-slack column
                %   is the transfer from its bus to the normalized weights, e_j - w, whose
                %   columns sum to zero, so makePTDF's txfr computes it column by column
                txfr = full(sparse(buses, 1:numel(buses), 1, nb, numel(buses))) ...
                    - (w / sum(w)) * ones(1, numel(buses));
                Hd = makePTDF(mpc, slack, txfr);
            end
            write_matrix(fullfile(out_dir, [name '_ptdf_distributed.csv']), labels(buses), ...
                rows, Hd);
        end
        if nb <= 300
            LODF = makeLODF(mpc, H);
        else
            % makeLODF's formula, column by column: the transfer from each outaged branch's
            %   from bus to its to bus, scaled by 1 / (1 - its own term), with -1 on the
            %   diagonal
            f = mpc.branch(outages, F_BUS);
            t = mpc.branch(outages, T_BUS);
            k = numel(outages);
            txfr = full(sparse([f; t], [1:k 1:k]', [ones(k, 1); -ones(k, 1)], nb, k));
            Ht = makePTDF(mpc, slack, txfr);
            own = Ht(sub2ind(size(Ht), outages, (1:k)'));
            LODF = Ht ./ (1 - own');
            LODF(sub2ind(size(LODF), outages, (1:k)')) = -1;
        end
        write_matrix(fullfile(out_dir, [name '_lodf.csv']), rows(outages), rows, LODF);
    end

    % The optimal power flows: DC with MIPS forced, then AC with MATPOWER's default, MIPS
    for run = {{'dcopf', true}, {'opf', false}}
        [prefix, dc] = deal(run{1}{1}, run{1}{2});
        if ~wanted(prefix)
            continue;
        end
        try
            if dc
                [result, success] = rundcopf(name, mpoption(mpopt, 'opf.dc.solver', 'MIPS'));
            else
                [result, success] = runopf(name, mpopt);
            end
        catch err
            fprintf(2, '%s (%s): %s\n', name, prefix, err.message);
            failed{end + 1} = [name ' (' prefix ')'];
            continue;
        end
        output = result.raw.output;
        if ~strcmp(output.alg, 'MIPS')
            fprintf(2, '%s (%s): solved by %s, not MIPS\n', name, prefix, output.alg);
            failed{end + 1} = [name ' (' prefix ')'];
            continue;
        end
        rows = (1:size(result.branch, 1))';
        if dc
            write_csv(fullfile(out_dir, [name '_dcopf_bus.csv']), 'bus,va_deg,lam_p', ...
                '%d,%.17g,%.17g', result.bus(:, [BUS_I, VA, LAM_P]));
            write_csv(fullfile(out_dir, [name '_dcopf_gen.csv']), ...
                'row,bus,status,pg_mw,mu_pmin,mu_pmax', '%d,%d,%d,%.17g,%.17g,%.17g', ...
                [(1:size(result.gen, 1))', result.gen(:, [GEN_BUS, GEN_STATUS, PG, MU_PMIN, MU_PMAX])]);
            write_csv(fullfile(out_dir, [name '_dcopf_branch.csv']), ...
                'row,from,to,status,pf_mw,mu_sf,mu_st,mu_angmin,mu_angmax', ...
                '%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g', ...
                [rows, result.branch(:, [F_BUS, T_BUS, BR_STATUS, PF, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX])]);
        else
            write_csv(fullfile(out_dir, [name '_opf_bus.csv']), ...
                'bus,vm_pu,va_deg,lam_p,lam_q,mu_vmin,mu_vmax', ...
                '%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g', ...
                result.bus(:, [BUS_I, VM, VA, LAM_P, LAM_Q, MU_VMIN, MU_VMAX]));
            write_csv(fullfile(out_dir, [name '_opf_gen.csv']), ...
                'row,bus,status,pg_mw,qg_mvar,mu_pmin,mu_pmax,mu_qmin,mu_qmax', ...
                '%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g', ...
                [(1:size(result.gen, 1))', result.gen(:, [GEN_BUS, GEN_STATUS, PG, QG, ...
                    MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX])]);
            write_csv(fullfile(out_dir, [name '_opf_branch.csv']), ...
                'row,from,to,status,pf_mw,qf_mvar,pt_mw,qt_mvar,mu_sf,mu_st,mu_angmin,mu_angmax', ...
                '%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g', ...
                [rows, result.branch(:, [F_BUS, T_BUS, BR_STATUS, PF, QF, PT, QT, MU_SF, ...
                    MU_ST, MU_ANGMIN, MU_ANGMAX])]);
        end
        hist = output.hist;
        write_csv(fullfile(out_dir, [name '_' prefix '_history.csv']), ...
            'iteration,feascond,gradcond,compcond,costcond,gamma,stepsize,obj,alphap,alphad', ...
            ['%d' repmat(',%.17g', 1, 9)], [(0:numel(hist) - 1)', [hist.feascond]', ...
                [hist.gradcond]', [hist.compcond]', [hist.costcond]', [hist.gamma]', ...
                [hist.stepsize]', [hist.obj]', [hist.alphap]', [hist.alphad]']);
        % The summary last, so that its presence marks a complete run
        fid = fopen(fullfile(out_dir, [name '_' prefix '_summary.csv']), 'w');
        if dc
            % The default solver's optimum, as a check: GLPK when the costs are linear
            try
                default = rundcopf(name, mpopt);
                [default_f, default_alg] = deal(default.f, default.raw.output.alg);
                if ~default.success
                    default_f = NaN;
                end
            catch err
                fprintf(2, '%s (dcopf, default solver): %s\n', name, err.message);
                [default_f, default_alg] = deal(NaN, 'none');
            end
            fprintf(fid, 'objective,iterations,success,default_objective,default_solver\n');
            fprintf(fid, '%.17g,%d,%d,%.17g,%s\n', result.f, output.iterations, success, ...
                default_f, default_alg);
        else
            fprintf(fid, 'objective,iterations,success\n');
            fprintf(fid, '%.17g,%d,%d\n', result.f, output.iterations, success);
        end
        fclose(fid);
        printf('%s: %s %s in %d iterations, objective %.17g\n', name, prefix, ...
            output.message, output.iterations, result.f);
    end

    if ~wanted('dcpf')
        continue;
    end

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
        if ~wanted(prefix)
            continue;
        end
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
