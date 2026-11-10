function mpc = three_bus
% The C interface's test case, read by both its Rust and its C tests: three buses, two
% generators and three branches, with costs of both models and bus names. Every column holds
% values distinct from its neighbours', so a field read from the wrong column shows.

mpc.version = '2';
mpc.baseMVA = 100;

%% bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1.02	0	230	1	1.1	0.9;
	2	2	50	10	0	5	1	1	-2	230	2	1.1	0.9;
	3	1	80	20	1.5	0	2	0.98	-4	115	3	1.05	0.95;
];

%% bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	90	5	100	-100	1.02	100	1	200	10	11	12	13	14	15	16	17	18	19	20	21;
	2	40	3	50	-50	1	110	1	80	0	1.1	1.2	1.3	1.4	1.5	1.6	1.7	1.8	1.9	2	2.1;
];

%% fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.01	0.1	0.02	250	260	270	0	0	1	-360	360;
	2	3	0.02	0.2	0.04	150	160	170	1.05	-3	1	-30	30;
	1	3	0.03	0.25	0	0	0	0	0	0	0	-360	360;
];

%% model	startup	shutdown	n	values...
mpc.gencost = [
	2	1500	100	3	0.11	5	150	0	0	0;
	1	0	0	3	0	0	40	2000	80	4400;
];

mpc.bus_name = {
	'Alpha';
	'Beta';
	'Gamma';
};
