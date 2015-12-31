function mpc = case2Tree
% BASED ON																					    
% Power flow data for Burak Kocuk, ect. @ GT 
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7056568

fprintf('Case2Tree\n');
																					
%% MATPOWER Case Format : Version 2																					
mpc.version = '2';																					
																					
%%-----  Power Flow Data  -----%%																					
%% system MVA base																					
mpc.baseMVA = 100;																					
																					
%% bus data																					
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin								
mpc.bus = [																					
	1	2	75	-84.7	0	0	1	nan     nan     0	1	1.1	0.9;
	2	2	105	22.8    0	0	1	nan     nan     0	1	1.1	0.9;
];																					
																					
%% generator data																					
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [																					
	1	nan	nan	300	-30     nan     nan     1	250	75	0	0	0	0	0	0	0	0	0	0	0;
	2	nan	nan	300	-30     nan     nan     1	300	70	0	0	0	0	0	0	0	0	0	0	0;
];																					
																					
%% branch data																					
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax								
mpc.branch = [																					
	1	2	0.01008	0.0504	nan	nan	0	0	0	0	1	-360	360;								
];																					
																					
%%-----  OPF Data  -----%%																					
%% generator cost data																					
%	1	startup	shutdown	n	x1	y1	...	xn	yn												
%	2	startup	shutdown	n	c(n-1)	...	c0														

mpc.gencost = [																					
	2	0	0	2	5	0;														
	2	0	0	2	1.2	0;																											
];	

																			
																					
																					
																					
																					
																					
																					
																					
																					
