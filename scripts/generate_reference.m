% Generates DC power flow reference results by running MATPOWER's rundcpf. Run through
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
% Values are printed with 17 significant digits, enough to round-trip a double exactly.

% Marks this file as a script, so Octave accepts the function definition that follows. Octave
%   needs a script's functions defined before they are called.
1;

function write_csv(path, header, row_format, data)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', header);
    fprintf(fid, [row_format '\n'], data');
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
end

% The run's provenance, so a regeneration with a different MATPOWER is visible in the diff
fid = fopen(fullfile(out_dir, 'VERSIONS'), 'w');
fprintf(fid, 'MATPOWER %s\nOctave %s\n', mpver, OCTAVE_VERSION);
fclose(fid);

if ~isempty(failed)
    error('rundcpf failed on %d of %d cases: %s', numel(failed), numel(cases), ...
        strjoin(failed, ', '));
end
