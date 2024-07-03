% Synchrobus with four VTs, by T. Van Acker

function mpc = test_spf
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    155.0   1       1.00    1.00;
    2       3       0.000   0.000   0.00    0.00    1       1.00    0.00    155.0   1       1.00    1.00;
    100     1       1.000   1.000   0.00    0.00    1       1.00	0.00    35.4    1       1.00    1.00;
    200     1       1.000   1.000   0.00    0.00    1       1.00    0.00    35.4    1       1.00    1.00;
    1500    1       1.000   1.000   0.00    0.00    1       1.00	0.00    10.0    1       1.00    1.00;
    1600    1       1.000   1.000   0.00    0.00    1       1.00    0.00    10.0    1       1.00    1.00;
    1300    1       0.000   0.000   0.00    0.00    1       1.00	0.00    0.69    1       1.00    1.00;
];

%% generator data
%	bus	    Pg	    Qg	    Qmax    Qmin    Vg	    mBase	status  Pmax    Pmin
mpc.gen = [
	1	    0.0		0.0	    9999.0	-9999.0	1.00000	100.0	1	    9999.0	0.0; % 150kV-Left-Inf
    2	    0.0		0.0	    9999.0	-9999.0	1.00000	100.0	1	    9999.0	0.0; % 150kV-Right-Inf
];

%% branch data
%	fbus	tbus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
% 150kV - 35.4kV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	1	    100	    0.0	    0.1536	0.0	    125	    125	    125	    0.0	    0.0	    1	    -30.0	30.0; % 17 - VT3 (Snom = 125e6 VA, usc = 0.192)
	1	    200	    0.0	    0.1608	0.0	    125	    125	    125	    0.0	    0.0	    1	    -30.0	30.0; % 18 - VT4 (Snom = 125e6 VA, usc = 0.201)
	2	    1500    0.0	    0.1379	0.0	    145	    145	    145	    0.0	    0.0	    1	    -30.0	30.0; % 24 - VT15 (Snom = 145e6 VA, usc = 0.202)
	%2	    1600    0.0	    0.1379	0.0	    145	    145	    145	    0.0	    0.0	    1	    -30.0	30.0; % 25 - VT16 (Snom = 145e6 VA, usc = 0.202)
% 35.4kV - North %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	100		1300	0.0	    0.1708	0.0	    70	    70	    70	    0.0	    0.0	    1	    -30.0	30.0; % 27 - 35kV-BASF-R100 SB - single coil = 0.1708, double coil = 0.0559
	200		1300	0.0	    0.1708	0.0	    70	    70	    70	    0.0	    0.0	    1	    -30.0	30.0; % 28 - 35kV-BASF-R200-1 SB
	1500	1300	0.0	    0.1708	0.0	    70	    70	    70	    0.0	    0.0	    1	    -30.0	30.0; % 30 - 35kV-BASF-R1500 SB
	1600	1300	0.0	    0.1708	0.0	    70	    70	    70	    0.0	    0.0	    1	    -30.0	30.0; % 31 - 35kV-BASF-R1600 SB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
];
