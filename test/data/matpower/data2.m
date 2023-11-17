function mpc = case5_acdc()
%case 5 nodes    Power flow data for 5 bus, 2 generator case.
%   Please see 'help caseformat' for details on the case file format.
%
%   case file can be used together with dc case files "case5_stagg_....m"
%
%   Network data from ...
%   G.W. Stagg, A.H. El-Abiad, "Computer methods in power system analysis",
%   McGraw-Hill, 1968.
%
%   MATPOWER case file data provided by Jef Beerten.

%% MATPOWER Case Format : Version 1
%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm      Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1       3       0	0	0   0   1       1.06	0	345     1       1.1     0.9;
	2       2       20	10	0   0   1       1       0	345     1       1.1     0.9;
	3       1       45	15	0   0   1       1       0	345     1       1.1     0.9;
	4       1       40	5	0   0   1       1       0	345     1       1.1     0.9;
	5       1       60	10	0   0   1       1       0	345     1       1.1     0.9;
];

%% generator data
%	bus	Pg      Qg	Qmax	Qmin	Vg	mBase       status	Pmax	Pmin	pc1 pc2 qlcmin qlcmax qc2min qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
mpc.gen = [
	1	0       0	500      -500    1.06	100       1       250     10 0 0 0 0 0 0 0 0 0 0 0;
  2	40      0	300      -300    1      100       1       300     10 0 0 0 0 0 0 0 0 0 0 0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle
%	status angmin angmax
mpc.branch = [
    1   2   0.02    0.06    0.06    100   100   100     0       0       1 -60 60;
    1   3   0.08    0.24    0.05    100   100   100     0       0       1 -60 60;
    2   3   0.06    0.18    0.04    100   100   100     0       0       1 -60 60;
    2   4   0.06    0.18    0.04    100   100   100     0       0       1 -60 60;
    2   5   0.04    0.12    0.03    100   100   100     0       0       1 -60 60;
    3   4   0.01    0.03    0.02    100   100   100     0       0       1 -60 60;
    4   5   0.08    0.24    0.05    100   100   100     0       0       1 -60 60;								
	2   3   0.06    0.18    0.04    100   100   100     0       0       1 -60 60;								
	3	5	0.06    0.18    0.04    100   100   100     0       0       1 -60 60;							
	2   5   0.04    0.12    0.03    100   100   100     0       0       1 -60 60;
];


%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0  1	0;
	2	0	0	3 0	 2	0;
];

% adds current ratings to branch matrix
%column_names%	c_rating_a
mpc.branch_currents = [
100;100;100;100;100;100;100;100;100;100;
];