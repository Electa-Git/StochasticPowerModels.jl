%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  %%%%%
%%%%    IEEE PES Power Grid Library - Optimal Power Flow - v21.07     %%%%%
%%%%          (https://github.com/power-grid-lib/pglib-opf)           %%%%%
%%%%               Benchmark Group - Typical Operations               %%%%%
%%%%                         29 - July - 2021                         %%%%%
%%%%                                                                  %%%%%
%%%%   The file has been extended for the purpose of stochastic OPF.  %%%%%
%%%%                                                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   CASE5  Power flow data for modified 5 bus, 5 gen case based on PJM 5-bus system
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from ...
%     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
%     Proceedings of the 2010 IEEE Power & Energy Society General Meeting
%
%   Created by Rui Bo in 2006, modified in 2010, 2014.
%
%   Copyright (c) 2010 by The Institute of Electrical and Electronics Engineers (IEEE)
%   Licensed under the Creative Commons Attribution 4.0
%   International license, http://creativecommons.org/licenses/by/4.0/
%
%   Contact M.E. Brennan (me.brennan@ieee.org) for inquries on further reuse of
%   this dataset.
%
%	Notes for stochastic extension:
%	- a stochastic germ ω comprised of two sources of uncertainty:
%		~ Beta distribution 𝓑(2.0,2.0) for load at bus 2
%		~ Normal distribution 𝓝(0.0,1.0) for load at bus 3
%	- stochastic parameters
%		~ μ = Pd,
%		~ σ = 0.10 * μ, and 
%		~ ε = 0.15 → λ = quantile(𝓝(0.0,1.0), 1.0 - 0.15) = 1.03643
function mpc = pglib_opf_case5_pjm
mpc.version = '2';
mpc.baseMVA = 100.0;

%% area data
%	area	refbus
mpc.areas = [
	1	 	4;
];

%% bus data
%	bus_i	type	Pd		Qd		Gs		Bs		area	Vm		Va		baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 	2	 	0.0	 	0.0	 	0.0	 	0.0	 	1	    1.00000	0.00000	230.0	1	    1.10000	0.90000;
	2	 	1	 	300.0	98.61	0.0	 	0.0	 	1	    1.00000	0.00000	230.0	1	    1.10000	0.90000;
	3	 	2	 	300.0	98.61	0.0	 	0.0	 	1	    1.00000	0.00000	230.0	1	    1.10000	0.90000;
	4	 	3	 	400.0	131.47	0.0	 	0.0	 	1	    1.00000	0.00000	230.0	1	    1.10000	0.90000;
	5	 	2	 	0.0	 	0.0	 	0.0	 	0.0	 	1	    1.00000	0.00000	230.0	1	    1.10000	0.90000;
];

%column_names%  dst_id	μ       σ       λvmin   λvmax	
mpc.bus_sdata = [
                0		0.0		0.0		1.03643	1.03643; % 1
                1		300.0	30.0	1.03643	1.03643; % 2
                2		300.0	30.0	1.03643	1.03643; % 3
                0		400.0	0.0		1.03643	1.03643; % 4
                0		0.0		0.0		1.03643	1.03643; % 5
]

%% generator data
%	bus		Pg		Qg		Qmax	Qmin	Vg		mBase	status	Pmax	Pmin
mpc.gen = [
	1	 	20.0	0.0	 	127.5	-127.5  1.0	 	100.0	1	 	150.0	0.0;
	3	 	260.0	0.0	 	390.0	-390.0	1.0	 	100.0	1	 	430.0	0.0;
	4	 	100.0	0.0	 	150.0	-150.0	1.0	 	100.0	1	 	990.0	0.0;
	5	 	300.0	0.0	 	450.0	-450.0	1.0	 	100.0	1	 	150.0	0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   2.000000	  14.000000	   0.000000;
	2	 0.0	 0.0	 3	   3.000000	  11.000000	   0.000000;
	2	 0.0	 0.0	 3	   4.000000	  14.000000	   0.000000;
	2	 0.0	 0.0	 3	   1.000000	  13.000000	   0.000000;
];

%% generator stochastic data 
%column_names%  λpmin   λpmax   λqmin   λqmax
mpc.gen_sdata = [
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
]

%% branch data
%	fbus	tbus	r		x		b		rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 	2	 	0.00281	0.0281	0.00000	400.0	400.0	400.0	0.0	 	0.0	 	1	 	-30.0	30.0;
	1	 	4	 	0.00304	0.0304	0.00000	426.0	426.0	426.0	0.0	 	0.0	 	1	 	-30.0	30.0;
	1	 	5	 	0.00064	0.0064	0.00000	426.0	426.0	426.0	0.0	 	0.0	 	1	 	-30.0	30.0;
	2	 	3	 	0.00108	0.0108	0.00000	426.0	426.0	426.0	0.0	 	0.0	 	1	 	-30.0	30.0;
	3	 	4	 	0.00297	0.0297	0.00000	426.0	426.0	426.0	0.0	 	0.0	 	1	 	-30.0	30.0;
	4	 	5	 	0.00297	0.0297	0.00000	240.0	240.0	240.0	0.0	 	0.0	 	1	 	-30.0	30.0;
];

%% branch stochastic data 
%column_names%  λcmax
mpc.branch_sdata = [
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
]

%% stochastic data
%column_names%  dst         pa      pb
mpc.sdata = [
                'Beta'   	2.0     2.0;
                'Normal'    0.0     0.0;
]

% INFO    : === Translation Options ===
% INFO    : Phase Angle Bound:           30.0 (deg.)
% INFO    : Line Capacity Model:         stat
% INFO    : Setting Flat Start
% INFO    : Line Capacity PAB:           15.0 (deg.)
% INFO    : 
% INFO    : === Generator Bounds Update Notes ===
% INFO    : 
% INFO    : === Base KV Replacement Notes ===
% INFO    : 
% INFO    : === Transformer Setting Replacement Notes ===
% INFO    : 
% INFO    : === Line Capacity Stat Model Notes ===
% INFO    : Updated Thermal Rating: on line 1-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 426
% INFO    : Updated Thermal Rating: on line 1-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 426
% INFO    : Updated Thermal Rating: on line 2-3 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 426
% INFO    : Updated Thermal Rating: on line 3-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 426
% INFO    : 
% INFO    : === Line Capacity Monotonicity Notes ===
% INFO    : 
% INFO    : === Voltage Setpoint Replacement Notes ===
% INFO    : Bus 1	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 2	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 3	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 4	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 5	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : 
% INFO    : === Generator Setpoint Replacement Notes ===
% INFO    : Gen at bus 1	: Pg=40.0, Qg=0.0 -> Pg=20.0, Qg=0.0
% INFO    : Gen at bus 1	: Vg=1.0 -> Vg=1.0
% INFO    : Gen at bus 1	: Pg=170.0, Qg=0.0 -> Pg=85.0, Qg=0.0
% INFO    : Gen at bus 1	: Vg=1.0 -> Vg=1.0
% INFO    : Gen at bus 3	: Pg=323.49, Qg=0.0 -> Pg=260.0, Qg=0.0
% INFO    : Gen at bus 3	: Vg=1.0 -> Vg=1.0
% INFO    : Gen at bus 4	: Pg=0.0, Qg=0.0 -> Pg=100.0, Qg=0.0
% INFO    : Gen at bus 4	: Vg=1.0 -> Vg=1.0
% INFO    : Gen at bus 5	: Pg=466.51, Qg=0.0 -> Pg=300.0, Qg=0.0
% INFO    : Gen at bus 5	: Vg=1.0 -> Vg=1.0
% INFO    : 
% INFO    : === Writing Matpower Case File Notes ===
