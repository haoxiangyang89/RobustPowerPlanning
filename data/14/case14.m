%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  %%%%%
%%%%      NICTA Energy System Test Case Archive (NESTA) - v0.7.0      %%%%%
%%%%              Optimal Power Flow - Typical Operation              %%%%%
%%%%                         05 - June - 2017                         %%%%%
%%%%                                                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Power flow data for IEEE 14 bus test case.
%  This data was converted from IEEE Common Data Format
%  (ieee14cdf.txt) on 20-Sep-2004 by cdf2matp, rev. 1.11
%
%  Converted from IEEE CDF file from:
%       http://www.ee.washington.edu/research/pstca/
%
%  CDF Header:
%  08/19/93 UW ARCHIVE           100.0  1962 W IEEE 14 Bus Test Case
%
function mpc = nesta_case14_ieee
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	    0.00000	 0.0	 1	    1.06000	    0.94000;
	2	 2	 21.7	 12.7	 0.0	 0.0	 1	    1.04017	   -4.25303	 0.0	 1	    1.06000	    0.94000;
	3	 2	 94.2	 19.0	 0.0	 0.0	 1	    1.00788	  -12.18881	 0.0	 1	    1.06000	    0.94000;
	4	 1	 47.8	 -3.9	 0.0	 0.0	 1	    1.01043	   -9.73410	 0.0	 1	    1.06000	    0.94000;
	5	 1	 7.6	 1.6	 0.0	 0.0	 1	    1.01324	   -8.24923	 0.0	 1	    1.06000	    0.94000;
	6	 2	 11.2	 7.5	 0.0	 0.0	 1	    1.06000	  -13.80994	 0.0	 1	    1.06000	    0.94000;
	7	 1	 0.0	 0.0	 0.0	 0.0	 1	    1.04402	  -12.83516	 0.0	 1	    1.06000	    0.94000;
	8	 2	 0.0	 0.0	 0.0	 0.0	 1	    1.06000	  -12.83516	 0.0	 1	    1.06000	    0.94000;
	9	 1	 29.5	 16.6	 0.0	 19.0	 1	    1.04096	  -14.45372	 0.0	 1	    1.06000	    0.94000;
	10	 1	 9.0	 5.8	 0.0	 0.0	 1	    1.03679	  -14.63018	 0.0	 1	    1.06000	    0.94000;
	11	 1	 3.5	 1.8	 0.0	 0.0	 1	    1.04473	  -14.35082	 0.0	 1	    1.06000	    0.94000;
	12	 1	 6.1	 1.6	 0.0	 0.0	 1	    1.04467	  -14.67869	 0.0	 1	    1.06000	    0.94000;
	13	 1	 13.5	 5.8	 0.0	 0.0	 1	    1.03948	  -14.75038	 0.0	 1	    1.06000	    0.94000;
	14	 1	 14.9	 5.0	 0.0	 0.0	 1	    1.02209	  -15.60843	 0.0	 1	    1.06000	    0.94000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 208.317	 0.0	 10.0	 0.0	 1.06	 100.0	 1	 362	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 50.0	 0.0	 0.0	 0.0	 0.0;
	2	 63.0	 26.964	 32.0	 -32.0	 1.04017	 100.0	 1	 63	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 20.0	 0.0	 0.0	 0.0	 0.0;
	3	 0.0	 29.63	 40.0	 0.0	 1.00788	 100.0	 1	 0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	6	 0.0	 13.23	 24.0	 -6.0	 1.06	 100.0	 1	 0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
	8	 0.0	 9.618	 24.0	 -6.0	 1.06	 100.0	 1	 0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   0.000000	   1.026133	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	   0.480811	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.000000	   0.000000	   0.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.01938	 0.05917	 0.0528	 472	 472	 472	 0.0	 0.0	 1	 -30.0	 30.0;
	1	 5	 0.05403	 0.22304	 0.0492	 128	 128	 128	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 3	 0.04699	 0.19797	 0.0438	 145	 145	 145	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 4	 0.05811	 0.17632	 0.034	 158	 158	 158	 0.0	 0.0	 1	 -30.0	 30.0;
	2	 5	 0.05695	 0.17388	 0.0346	 161	 161	 161	 0.0	 0.0	 1	 -30.0	 30.0;
	3	 4	 0.06701	 0.17103	 0.0128	 160	 160	 160	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 5	 0.01335	 0.04211	 0.0	 664	 664	 664	 0.0	 0.0	 1	 -30.0	 30.0;
	4	 7	 0.0	 0.20912	 0.0	 141	 141	 141	 0.978	 0.0	 1	 -30.0	 30.0;
	4	 9	 0.0	 0.55618	 0.0	 53	 53	 53	 0.969	 0.0	 1	 -30.0	 30.0;
	5	 6	 0.0	 0.25202	 0.0	 117	 117	 117	 0.932	 0.0	 1	 -30.0	 30.0;
	6	 11	 0.09498	 0.1989	 0.0	 134	 134	 134	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 12	 0.12291	 0.25581	 0.0	 104	 104	 104	 0.0	 0.0	 1	 -30.0	 30.0;
	6	 13	 0.06615	 0.13027	 0.0	 201	 201	 201	 0.0	 0.0	 1	 -30.0	 30.0;
	7	 8	 0.0	 0.17615	 0.0	 167	 167	 167	 0.0	 0.0	 1	 -30.0	 30.0;
	7	 9	 0.0	 0.11001	 0.0	 267	 267	 267	 0.0	 0.0	 1	 -30.0	 30.0;
	9	 10	 0.03181	 0.0845	 0.0	 325	 325	 325	 0.0	 0.0	 1	 -30.0	 30.0;
	9	 14	 0.12711	 0.27038	 0.0	 99	 99	 99	 0.0	 0.0	 1	 -30.0	 30.0;
	10	 11	 0.08205	 0.19207	 0.0	 141	 141	 141	 0.0	 0.0	 1	 -30.0	 30.0;
	12	 13	 0.22092	 0.19988	 0.0	 99	 99	 99	 0.0	 0.0	 1	 -30.0	 30.0;
	13	 14	 0.17093	 0.34802	 0.0	 76	 76	 76	 0.0	 0.0	 1	 -30.0	 30.0;
];

% INFO    : === Translation Options ===
% INFO    : Phase Angle Bound:           30.0 (deg.)
% INFO    : Line Capacity Model:         stat
% INFO    : Gen Active Capacity Model:   stat
% INFO    : Gen Reactive Capacity Model: am50ag
% INFO    : Gen Active Cost Model:       stat
% INFO    : AC OPF Solution File:        nesta_case14_ieee.m.opf.sol
% INFO    : Line Capacity PAB:           15.0 (deg.)
% INFO    :
% INFO    : === Generator Classification Notes ===
% INFO    : SYNC   3   -     0.00
% INFO    : NG     2   -   100.00
% INFO    :
% INFO    : === Generator Active Capacity Stat Model Notes ===
% INFO    : Gen at bus 1 - NG	: Pg=232.4, Pmax=332.4 -> Pmax=362   samples: 15
% INFO    : Gen at bus 2 - NG	: Pg=40.0, Pmax=140.0 -> Pmax=63   samples: 1
% INFO    : Gen at bus 3 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    : Gen at bus 6 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    : Gen at bus 8 - SYNC	: Pg=0.0, Pmax=100.0 -> Pmax=0   samples: 0
% INFO    :
% INFO    : === Generator Reactive Capacity Atmost Max 50 Percent Active Model Notes ===
% INFO    : Gen at bus 2 - NG	: Pmax 63.0, Qmin -40.0, Qmax 50.0 -> Qmin -32.0, Qmax 32.0
% INFO    :
% INFO    : === Generator Active Cost Stat Model Notes ===
% INFO    : Updated Generator Cost: NG - 0.0 20.0 0.0430293 -> 0 1.02613295282 0
% INFO    : Updated Generator Cost: NG - 0.0 20.0 0.25 -> 0 0.480811431055 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    : Updated Generator Cost: SYNC - 0.0 40.0 0.01 -> 0 0.0 0
% INFO    :
% INFO    : === Line Capacity Stat Model Notes ===
% WARNING : Missing data for branch flow stat model on line 1-2 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.01938 x=0.05917
% INFO    : Updated Thermal Rating: on line 1-2 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 472
% WARNING : Missing data for branch flow stat model on line 1-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05403 x=0.22304
% INFO    : Updated Thermal Rating: on line 1-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 128
% WARNING : Missing data for branch flow stat model on line 2-3 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.04699 x=0.19797
% INFO    : Updated Thermal Rating: on line 2-3 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 145
% WARNING : Missing data for branch flow stat model on line 2-4 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05811 x=0.17632
% INFO    : Updated Thermal Rating: on line 2-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 158
% WARNING : Missing data for branch flow stat model on line 2-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.05695 x=0.17388
% INFO    : Updated Thermal Rating: on line 2-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 161
% WARNING : Missing data for branch flow stat model on line 3-4 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.06701 x=0.17103
% INFO    : Updated Thermal Rating: on line 3-4 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 160
% WARNING : Missing data for branch flow stat model on line 4-5 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.01335 x=0.04211
% INFO    : Updated Thermal Rating: on line 4-5 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 664
% WARNING : Missing data for branch flow stat model on line 4-7 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.20912
% INFO    : Updated Thermal Rating: on transformer 4-7 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 141
% WARNING : Missing data for branch flow stat model on line 4-9 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.55618
% INFO    : Updated Thermal Rating: on transformer 4-9 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 53
% WARNING : Missing data for branch flow stat model on line 5-6 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.25202
% INFO    : Updated Thermal Rating: on transformer 5-6 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 117
% WARNING : Missing data for branch flow stat model on line 6-11 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.09498 x=0.1989
% INFO    : Updated Thermal Rating: on line 6-11 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 134
% WARNING : Missing data for branch flow stat model on line 6-12 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.12291 x=0.25581
% INFO    : Updated Thermal Rating: on line 6-12 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 104
% WARNING : Missing data for branch flow stat model on line 6-13 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.06615 x=0.13027
% INFO    : Updated Thermal Rating: on line 6-13 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 201
% WARNING : Missing data for branch flow stat model on line 7-8 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.17615
% INFO    : Updated Thermal Rating: on line 7-8 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 167
% WARNING : Missing data for branch flow stat model on line 7-9 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.0 x=0.11001
% INFO    : Updated Thermal Rating: on line 7-9 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 267
% WARNING : Missing data for branch flow stat model on line 9-10 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.03181 x=0.0845
% INFO    : Updated Thermal Rating: on line 9-10 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 325
% WARNING : Missing data for branch flow stat model on line 9-14 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.12711 x=0.27038
% INFO    : Updated Thermal Rating: on line 9-14 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 99
% WARNING : Missing data for branch flow stat model on line 10-11 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.08205 x=0.19207
% INFO    : Updated Thermal Rating: on line 10-11 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 141
% WARNING : Missing data for branch flow stat model on line 12-13 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.22092 x=0.19988
% INFO    : Updated Thermal Rating: on line 12-13 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 99
% WARNING : Missing data for branch flow stat model on line 13-14 using max current model : from_basekv=0.0 to_basekv=0.0 r=0.17093 x=0.34802
% INFO    : Updated Thermal Rating: on line 13-14 : Rate A, Rate B, Rate C , 9900.0, 0.0, 0.0 -> 76
% INFO    :
% INFO    : === Voltage Setpoint Replacement Notes ===
% INFO    : Bus 1	: V=1.06, theta=0.0 -> V=1.06, theta=0.0
% INFO    : Bus 2	: V=1.045, theta=-4.98 -> V=1.04017, theta=-4.25303
% INFO    : Bus 3	: V=1.01, theta=-12.72 -> V=1.00788, theta=-12.18881
% INFO    : Bus 4	: V=1.019, theta=-10.33 -> V=1.01043, theta=-9.7341
% INFO    : Bus 5	: V=1.02, theta=-8.78 -> V=1.01324, theta=-8.24923
% INFO    : Bus 6	: V=1.07, theta=-14.22 -> V=1.06, theta=-13.80994
% INFO    : Bus 7	: V=1.062, theta=-13.37 -> V=1.04402, theta=-12.83516
% INFO    : Bus 8	: V=1.09, theta=-13.36 -> V=1.06, theta=-12.83516
% INFO    : Bus 9	: V=1.056, theta=-14.94 -> V=1.04096, theta=-14.45372
% INFO    : Bus 10	: V=1.051, theta=-15.1 -> V=1.03679, theta=-14.63018
% INFO    : Bus 11	: V=1.057, theta=-14.79 -> V=1.04473, theta=-14.35082
% INFO    : Bus 12	: V=1.055, theta=-15.07 -> V=1.04467, theta=-14.67869
% INFO    : Bus 13	: V=1.05, theta=-15.16 -> V=1.03948, theta=-14.75038
% INFO    : Bus 14	: V=1.036, theta=-16.04 -> V=1.02209, theta=-15.60843
% INFO    :
% INFO    : === Generator Setpoint Replacement Notes ===
% INFO    : Gen at bus 1	: Pg=232.4, Qg=-16.9 -> Pg=208.317, Qg=0.0
% INFO    : Gen at bus 1	: Vg=1.06 -> Vg=1.06
% INFO    : Gen at bus 2	: Pg=40.0, Qg=42.4 -> Pg=63.0, Qg=26.964
% INFO    : Gen at bus 2	: Vg=1.045 -> Vg=1.04017
% INFO    : Gen at bus 3	: Pg=0.0, Qg=23.4 -> Pg=0.0, Qg=29.63
% INFO    : Gen at bus 3	: Vg=1.01 -> Vg=1.00788
% INFO    : Gen at bus 6	: Pg=0.0, Qg=12.2 -> Pg=0.0, Qg=13.23
% INFO    : Gen at bus 6	: Vg=1.07 -> Vg=1.06
% INFO    : Gen at bus 8	: Pg=0.0, Qg=17.4 -> Pg=0.0, Qg=9.618
% INFO    : Gen at bus 8	: Vg=1.09 -> Vg=1.06
% INFO    :
% INFO    : === Writing Matpower Case File Notes ===
