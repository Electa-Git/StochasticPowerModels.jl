%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  %%%%%
%%%%    IEEE PES Power Grid Library - Optimal Power Flow - v19.01     %%%%%
%%%%          (https://github.com/power-grid-lib/pglib-opf)           %%%%%
%%%%               Benchmark Group - Typical Operations               %%%%%
%%%%                         04 - Sep - 2021                          %%%%%
%%%%                                                                  %%%%%
%%%%   The file has been extended for the purpose of stochastic OPF.  %%%%%
%%%%                                                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Power flow data for 30 bus, 6 generator case.
%
%   Based on data from:
%     Alsac, O. & Stott, B., "Optimal Load Flow with Steady State Security",
%     IEEE Transactions on Power Apparatus and Systems, Vol. PAS 93, No. 3,
%     1974, pp. 745-751.
%   
%   The generation MVA limit is used for Q upper bounds
%   Shunts assumed to be in real value (not p.u. as specified in the paper)
%   
%   Additional modifications based on:
%     Ferrero, R.W., Shahidehpour, S.M., Ramesh, V.C., "Transaction analysis
%     in deregulated power systems using game theory", IEEE Transactions on
%     Power Systems, Vol. 12, No. 3, Aug 1997, pp. 1340-1347.
%
%   Generators moved to match Figure 3
%   Removed load at bus 5 (Figure 3)
%   Costs and Generator bounds updated to Table 1
%
%   Copyright (c) 1997 by The Institute of Electrical and Electronics Engineers (IEEE)
%   Licensed under the Creative Commons Attribution 4.0
%   International license, http://creativecommons.org/licenses/by/4.0/
%
%   Contact M.E. Brennan (me.brennan@ieee.org) for inquries on further reuse of
%   this dataset.
%
%   Notes for stochastic extension:
%	- shunt elements are neglected
% 	- voltage magnitude at slack bus 1 is assumed constant at one
%	- line current limits are set to the nominal values of the per-unit line 
%	  ratings, except for two lines where the capacity is reduced from 16 to, 
%	  15-23: imax = 16 -> 11, and 25-27: imax = 16 -> 12
% 	- a stochastic germ ω comprised of four distinct sources of uncertainty, two
%	  Beta distributions (one symmetric, one non-symmetric) and two normal 
%	  distributions as described in Table VI. [1]
% 	- additional stochastic parameters: 
%		- σ = 0.15 * μ, and 
%		- ε = 0.15: quantile(𝓝(0.0,1.0), 1.0 - 0.15) = 1.03643
% 	
% 	[1] Chance-constrained AC optimal power flow - a polynomial chaos approach 
%		by Muhlpfordt et al.


function mpc = pglib_opf_case30_fsr
mpc.version = '2';
mpc.baseMVA = 100.0;

%% area data
%	area	refbus
mpc.areas = [
	1	 8;
	2	 23;
	3	 26;
];

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	135	1	1.1	0.95;
	2	2	21.7	12.7	0	0	1	1	0	135	1	1.1	0.95;
	3	1	2.4	1.2	0	0	1	1	0	135	1	1.05	0.95;
	4	1	7.6	1.6	0	0	1	1	0	135	1	1.05	0.95;
	5	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	6	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	7	1	22.8	10.9	0	0	1	1	0	135	1	1.05	0.95;
	8	1	30	30	0	0	1	1	0	135	1	1.05	0.95;
	9	1	0.00	0.00	0	0	1	1	0	135	1	1.05	0.95;
	10	1	5.80	2.00	0	0	3	1	0	135	1	1.05	0.95;
	11	1	6.00	2.20	0	0	1	1	0	135	1	1.05	0.95;
	12	1	11.2	7.5	0	0	2	1	0	135	1	1.05	0.95;
	13	2	0	0	0	0	2	1	0	135	1	1.1	0.95;
	14	1	6.2	1.6	0	0	2	1	0	135	1	1.05	0.95;
	15	1	8.2	2.5	0	0	2	1	0	135	1	1.05	0.95;
	16	1	3.5	1.8	0	0	2	1	0	135	1	1.05	0.95;
	17	1	9	5.8	0	0	2	1	0	135	1	1.05	0.95;
	18	1	3.2	0.9	0	0	2	1	0	135	1	1.05	0.95;
	19	1	9.5	3.4	0	0	2	1	0	135	1	1.05	0.95;
	20	1	2.2	0.7	0	0	2	1	0	135	1	1.05	0.95;
	21	1	17.5	11.2	0	0	3	1	0	135	1	1.05	0.95;
	22	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	23	2	3.2	1.6	0	0	2	1	0	135	1	1.1	0.95;
	24	1	8.7	6.7	0	0	3	1	0	135	1	1.05	0.95;
	25	1	0	0	0	0	3	1	0	135	1	1.05	0.95;
	26	1	3.5	2.3	0	0	3	1	0	135	1	1.05	0.95;
	27	2	0	0	0	0	3	1	0	135	1	1.1	0.95;
	28	1	0	0	0	0	1	1	0	135	1	1.05	0.95;
	29	1	2.4	0.9	0	0	3	1	0	135	1	1.05	0.95;
	30	1	10.6	1.9	0	0	3	1	0	135	1	1.05	0.95;
];

%column_names%  dst_id  μ       σ       λvmin   λvmax
mpc.bus_sdata = [
                0       0.0     0.0     1.03643 1.03643; % 1
                1       21.7    3.255   1.03643 1.03643; % 2
                1       2.4     0.36    1.03643 1.03643; % 3
                2       7.6     1.14    1.03643 1.03643; % 4
                0       0.0    	0.0    	1.03643 1.03643; % 5
                0       0.0     0.0     1.03643 1.03643; % 6
                0       22.8    0.0     1.03643 1.03643; % 7
                0       30.0    0.0     1.03643 1.03643; % 8
                0       0.0     0.0     1.03643 1.03643; % 9 
                4       5.8     0.87    1.03643 1.03643; % 10
                0       6.0     0.0     1.03643 1.03643; % 11
                0       11.2    0.0     1.03643 1.03643; % 12
                0       0.0     0.0     1.03643 1.03643; % 13
                0       6.2     0.0     1.03643 1.03643; % 14 
                0       8.2     0.0     1.03643 1.03643; % 15
                0       3.5     0.0     1.03643 1.03643; % 16
                0       9.0     0.0     1.03643 1.03643; % 17
                0       3.2     0.0     1.03643 1.03643; % 18
                0       9.5     0.0     1.03643 1.03643; % 19
                0       2.2     0.0     1.03643 1.03643; % 20
                4       17.5    2.625   1.03643 1.03643; % 21
                0       0.0     0.0     1.03643 1.03643; % 22
                0       3.2     0.0     1.03643 1.03643; % 23
                3       8.7     1.305   1.03643 1.03643; % 24
                0       0.0     0.0     1.03643 1.03643; % 25
                0       3.5     0.0     1.03643 1.03643; % 26
                0       0.0     0.0     1.03643 1.03643; % 27
                0       0.0     0.0     1.03643 1.03643; % 28
                0       2.4     0.0     1.03643 1.03643; % 29
                0       10.6    0.0     1.03643 1.03643; % 30
]

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	1	23.54	0	150.0	-20	1	100	1	80	0;
	2	60.97	0	60.0	-20	1	100	1	58	0;
	22	21.59	0	62.5	-15	1	100	1	24	0;
	27	26.91	0	48.7	-15	1	100	1	35	0;
	23	19.20	0	40.0	-10	1	100	1	30	0;
	13	37.00	0	44.7	-15	1	100	1	40	0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	0.0200	2.00	0;
	2	0	0	3	0.0175	1.75	0;
	2	0	0	3	0.0625	1.00	0;
	2	0	0	3	0.00834	3.25	0;
	2	0	0	3	0.0250	3.00	0;
	2	0	0	3	0.0250	3.00	0;
];

%% generator stochastic data 
%column_names%  λpmin   λpmax   λqmin   λqmax
mpc.gen_sdata = [
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
                1.03643 1.03643 1.03643 1.03643;
]

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.02	0.06	0.00	130	130	130	0	0	1	-360	360;
	1	3	0.05	0.19	0.00	130	130	130	0	0	1	-360	360;
	2	4	0.06	0.17	0.00	65	65	65	0	0	1	-360	360;
	3	4	0.01	0.04	0.00	130	130	130	0	0	1	-360	360;
	2	5	0.05	0.20	0.00	130	130	130	0	0	1	-360	360;
	2	6	0.06	0.18	0.00	65	65	65	0	0	1	-360	360;
	4	6	0.01	0.04	0.00	90	90	90	0	0	1	-360	360;
	5	7	0.05	0.12	0.00	70	70	70	0	0	1	-360	360;
	6	7	0.03	0.08	0.00	130	130	130	0	0	1	-360	360;
	6	8	0.01	0.04	0.00	35	32	32	0	0	1	-360	360;
	6	9	0.00	0.21	0.00	65	65	65	0	0	1	-360	360;
	6	10	0.00	0.56	0.00	32	32	32	0	0	1	-360	360;
	9	11	0.00	0.21	0.00	65	65	65	0	0	1	-360	360;
	9	10	0.00	0.11	0.00	65	65	65	0	0	1	-360	360;
	4	12	0.00	0.26	0.00	65	65	65	0	0	1	-360	360;
	12	13	0.00	0.14	0.00	65	65	65	0	0	1	-360	360;
	12	14	0.12	0.26	0.00	32	32	32	0	0	1	-360	360;
	12	15	0.07	0.13	0.00	32	32	32	0	0	1	-360	360;
	12	16	0.09	0.20	0.00	32	32	32	0	0	1	-360	360;
	14	15	0.22	0.20	0.00	14	16	16	0	0	1	-360	360;
	16	17	0.08	0.19	0.00	16	16	16	0	0	1	-360	360;
	15	18	0.11	0.22	0.00	16	16	16	0	0	1	-360	360;
	18	19	0.06	0.13	0.00	16	16	16	0	0	1	-360	360;
	19	20	0.03	0.07	0.00	32	32	32	0	0	1	-360	360;
	10	20	0.09	0.21	0.00	32	32	32	0	0	1	-360	360;
	10	17	0.03	0.08	0.00	32	32	32	0	0	1	-360	360;
	10	21	0.03	0.07	0.00	32	32	32	0	0	1	-360	360;
	10	22	0.07	0.15	0.00	32	32	32	0	0	1	-360	360;
	21	22	0.01	0.02	0.00	32	32	32	0	0	1	-360	360;
	15	23	0.10	0.20	0.00	11	16	16	0	0	1	-360	360;
	22	24	0.12	0.18	0.00	16	16	16	0	0	1	-360	360;
	23	24	0.13	0.27	0.00	16	16	16	0	0	1	-360	360;
	24	25	0.19	0.33	0.00	16	16	16	0	0	1	-360	360;
	25	26	0.25	0.38	0.00	16	16	16	0	0	1	-360	360;
	25	27	0.11	0.21	0.00	12	16	16	0	0	1	-360	360;
	28	27	0.00	0.40	0.00	65	65	65	0	0	1	-360	360;
	27	29	0.22	0.42	0.00	16	16	16	0	0	1	-360	360;
	27	30	0.32	0.60	0.00	16	16	16	0	0	1	-360	360;
	29	30	0.24	0.45	0.00	16	16	16	0	0	1	-360	360;
	8	28	0.06	0.20	0.00	32	32	32	0	0	1	-360	360;
	6	28	0.02	0.06	0.00	32	32	32	0	0	1	-360	360;
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
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
                1.03643;
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
                'Beta'      2.0     5.0;
                'Normal'    0.0     0.0;
				'Normal'	0.0		0.0;
]

% INFO    : === Translation Options ===
% INFO    : Phase Angle Bound:           30.0 (deg.)
% INFO    : Setting Flat Start
% INFO    : 
% INFO    : === Base KV Replacement Notes ===
% INFO    : 
% INFO    : === Transformer Setting Replacement Notes ===
% INFO    : 
% INFO    : === Line Capacity Monotonicity Notes ===
% INFO    : 
% INFO    : === Voltage Setpoint Replacement Notes ===
% INFO    : Bus 1	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 2	: V=1.0, theta=0.0 -> V=1.025, theta=0.0
% INFO    : Bus 3	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 4	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 5	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 6	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 7	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 8	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 9	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 10	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 11	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 12	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 13	: V=1.0, theta=0.0 -> V=1.025, theta=0.0
% INFO    : Bus 14	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 15	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 16	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 17	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 18	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 19	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 20	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 21	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 22	: V=1.0, theta=0.0 -> V=1.025, theta=0.0
% INFO    : Bus 23	: V=1.0, theta=0.0 -> V=1.025, theta=0.0
% INFO    : Bus 24	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 25	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 26	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 27	: V=1.0, theta=0.0 -> V=1.025, theta=0.0
% INFO    : Bus 28	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 29	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : Bus 30	: V=1.0, theta=0.0 -> V=1.0, theta=0.0
% INFO    : 
% INFO    : === Generator Setpoint Replacement Notes ===
% INFO    : Gen at bus 1	: Pg=23.54, Qg=0.0 -> Pg=40.0, Qg=115.0
% INFO    : Gen at bus 1	: Vg=1.0 -> Vg=1.0
% INFO    : Gen at bus 2	: Pg=60.97, Qg=0.0 -> Pg=40.0, Qg=40.0
% INFO    : Gen at bus 2	: Vg=1.0 -> Vg=1.025
% INFO    : Gen at bus 22	: Pg=21.59, Qg=0.0 -> Pg=25.0, Qg=32.5
% INFO    : Gen at bus 22	: Vg=1.0 -> Vg=1.025
% INFO    : Gen at bus 27	: Pg=26.91, Qg=0.0 -> Pg=27.5, Qg=22.5
% INFO    : Gen at bus 27	: Vg=1.0 -> Vg=1.025
% INFO    : Gen at bus 23	: Pg=19.2, Qg=0.0 -> Pg=15.0, Qg=20.0
% INFO    : Gen at bus 23	: Vg=1.0 -> Vg=1.025
% INFO    : Gen at bus 13	: Pg=37.0, Qg=0.0 -> Pg=20.0, Qg=22.5
% INFO    : Gen at bus 13	: Vg=1.0 -> Vg=1.025
% INFO    : 
% INFO    : === Writing Matpower Case File Notes ===
