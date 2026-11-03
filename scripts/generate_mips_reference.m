% Records MIPS's results on the test problems of MATPOWER's own MIPS tests (mips/lib/t/t_mips.m
%   and t_qps_mips.m), restated here, for interior_point's tests: each problem solved with step
%   control off and on, printing the exit flag, iteration count, objective, solution and
%   multipliers with 17 significant digits. Run through scripts/generate-mips-reference.sh, which
%   fetches MATPOWER and writes interior_point/tests/data/mips_reference.txt.
%
% Each line is: name, step_control, then key=values pairs separated by '|', with vectors as
%   space-separated numbers (empty when MIPS returns none).
%
% Usage: octave generate_mips_reference.m <matpower-dir>

1;

function out = vec(v)
    out = strtrim(sprintf('%.17g ', v));
end

function report(name, sc, x, f, eflag, output, lam)
    fields = {};
    fields{end + 1} = sprintf('eflag=%d', eflag);
    fields{end + 1} = sprintf('iterations=%d', output.iterations);
    fields{end + 1} = sprintf('f=%.17g', f);
    fields{end + 1} = ['x=' vec(x)];
    for name_ = {'eqnonlin', 'ineqnonlin', 'mu_l', 'mu_u', 'lower', 'upper'}
        if isfield(lam, name_{1})
            fields{end + 1} = [name_{1} '=' vec(lam.(name_{1}))];
        end
    end
    printf('%s,%d,%s\n', name, sc, strjoin(fields, '|'));
end

% The problems of t_mips.m: objective, constraints and Hessians
function [f, df, d2f] = f2(x)
    a = 100;
    f = a*(x(2)-x(1)^2)^2+(1-x(1))^2;
    df = [  4*a*(x(1)^3 - x(1)*x(2)) + 2*x(1)-2;
            2*a*(x(2) - x(1)^2)                     ];
    d2f = 4*a*[ 3*x(1)^2 - x(2) + 1/(2*a),  -x(1);
                -x(1)                       1/2       ];
end
function [f, df, d2f] = f3(x)
    H = [5 -2 -1; -2 4 3; -1 3 5];
    c = [2; -35; -47];
    f = 1/2 * x'*H*x + c'*x + 5;
    df = H*x + c;
    d2f = H;
end
function [f, df, d2f] = f4(x)
    H = [1003.1 4.3 6.3 5.9; 4.3 2.2 2.1 3.9; 6.3 2.1 3.5 4.8; 5.9 3.9 4.8 10];
    f = 1/2 * x'*H*x;
    df = H*x;
    d2f = H;
end
function [f, df, d2f] = f5(x)
    c = -[1; 1];
    f = c'*x;
    df = c;
    d2f = zeros(2,2);
end
function [h, g, dh, dg] = gh5(x)
    h = [ -1 -1; 1 1] * x.^2 + [1; -2];
    dh = 2 * [-x(1) x(1); -x(2) x(2)];
    g = []; dg = [];
end
function Lxx = hess5(x, lam, cost_mult)
    mu = lam.ineqnonlin;
    Lxx = 2*[-1 1]*mu*eye(2);
end
function [f, df, d2f] = f6(x)
    f = -x(1)*x(2) - x(2)*x(3);
    df = -[x(2); x(1)+x(3); x(2)];
    d2f = -[0 1 0; 1 0 1; 0 1 0];
end
function [h, g, dh, dg] = gh6(x)
    h = [ 1 -1 1; 1 1 1] * x.^2 + [-2; -10];
    dh = 2 * [x(1) x(1); -x(2) x(2); x(3) x(3)];
    g = []; dg = [];
end
function Lxx = hess6(x, lam, cost_mult)
    mu = lam.ineqnonlin;
    Lxx = cost_mult * [0 -1 0; -1 0 -1; 0 -1 0] + ...
            [2*[1 1]*mu 0 0; 0 2*[-1 1]*mu 0; 0 0 2*[1 1]*mu];
end
function [f, df, d2f] = f7(x)
    f = x(1)*x(4)*sum(x(1:3)) + x(3);
    df = [ x(1)*x(4) + x(4)*sum(x(1:3)); x(1)*x(4); x(1)*x(4) + 1; x(1)*sum(x(1:3)) ];
    d2f = sparse([ 2*x(4) x(4) x(4) 2*x(1)+x(2)+x(3); x(4) 0 0 x(1); x(4) 0 0 x(1); ...
                   2*x(1)+x(2)+x(3) x(1) x(1) 0 ]);
end
function [h, g, dh, dg] = gh7(x)
    g = sum(x.^2) - 40;
    h = -prod(x) + 25;
    dg = 2*x;
    dh = -prod(x)./x;
end
function Lxx = hess7(x, lam, sigma)
    lambda = lam.eqnonlin;
    mu     = lam.ineqnonlin;
    [f, df, d2f] = f7(x);
    Lxx = sigma * d2f + lambda*2*speye(4) - ...
       mu*sparse([ 0 x(3)*x(4) x(2)*x(4) x(2)*x(3); x(3)*x(4) 0 x(1)*x(4) x(1)*x(3); ...
                   x(2)*x(4) x(1)*x(4) 0 x(1)*x(2); x(2)*x(3) x(1)*x(3) x(1)*x(2) 0 ]);
end

args = argv();
addpath(genpath(args{1}));
for sc = [0 1]
    opt = struct('verbose', 0, 'step_control', sc);
    [x, f, s, out, lam] = mips(@f2, [-1.9; 2], [], [], [], [], [], [], [], opt);
    report('banana', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = mips(@f3, [0; 0; 0], [], [], [], [], [], [], [], opt);
    report('quadratic3', sc, x, f, s, out, lam);
    A = [1 1 1 1; 0.17 0.11 0.10 0.18];
    [x, f, s, out, lam] = mips(@f4, [1; 0; 0; 1], A, [1; 0.10], [1; Inf], zeros(4,1), [], [], [], opt);
    report('qp4', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = mips(@f5, [1.1; 0], [], [], [], zeros(2,1), [], @gh5, @hess5, opt);
    report('nonlinear2', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = mips(@f6, [1; 1; 0], [], [], [], [], [], @gh6, @hess6, opt);
    report('nonlinear3', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = mips(@f7, [1; 5; 5; 1], [], [], [], ones(4,1), 5*ones(4,1), @gh7, @hess7, opt);
    report('hs71', sc, x, f, s, out, lam);

    % The problems of t_qps_mips.m, which sets comptol 1e-8
    qopt = opt;
    qopt.comptol = 1e-8;
    [x, f, s, out, lam] = qps_mips([], [-5; -4; -6], [1 -1 1; -3 -2 -4; 3 2 0], ...
        [-Inf; -42; -Inf], [20; Inf; 30], [0; 0; 0], [], [], qopt);
    report('lp3', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = qps_mips([5 -2 -1; -2 4 3; -1 3 5], [2; -35; -47], [], [], [], ...
        [], [], [0; 0; 0], qopt);
    report('qps_quadratic3', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = qps_mips([1 -1; -1 2], [-2; -6], [1 1; -1 2; 2 1], [], [2; 2; 3], ...
        [0; 0], [], [], qopt);
    report('qp2', sc, x, f, s, out, lam);
    [x, f, s, out, lam] = qps_mips([1003.1 4.3 6.3 5.9; 4.3 2.2 2.1 3.9; 6.3 2.1 3.5 4.8; ...
        5.9 3.9 4.8 10], zeros(4,1), [1 1 1 1; 0.17 0.11 0.10 0.18], [1; 0.10], [1; Inf], ...
        zeros(4,1), [], [1; 0; 0; 1], qopt);
    report('qps_qp4', sc, x, f, s, out, lam);
    p = struct('A', sparse([1 1]), 'c', [1;1], 'u', -1, 'xmin', [0;0], 'opt', qopt);
    [x, f, s, out, lam] = qps_mips(p);
    report('infeasible_lp', sc, x, f, s, out, lam);
end
