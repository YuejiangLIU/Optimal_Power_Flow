function mpc = case3Tree
% BASED ON																					    
% Power flow data for Burak Kocuk, ect. @ GT 
% http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7056568

fprintf('Case3Tree\n');

%% MATPOWER Case Format : Version 2																					
mpc.version = '2';																					
																					
%%-----  Power Flow Data  -----%%																					
%% system MVA base																					
mpc.baseMVA = 100;																					
																					
%% bus data																					
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin								
mpc.bus = [																					
	1	2	50	-52.3	0	0	1	nan     nan     0	1	1.1	0.9;
	2	1	70	14.1	0	0	1	nan     nan     0	1	1.1	0.9;
	3	1	60	-82.3	0	0	1	nan     nan     0	1	1.1	0.9;
];																					
																					
%% generator data																					
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [																					
	1	nan	nan	500 	-100    nan nan     1	550     150     0	0	0	0	0	0	0	0	0	0	0;
];																					
																					
%% branch data																					
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax								
mpc.branch = [																					
	1	2	0.01008	0.0504	nan	nan	0	0	0	0	1	-360	360;								
	2	3	0.07500	0.0840	nan	nan	0	0	0	0	1	-360	360;
];																					
																					
%%-----  OPF Data  -----%%																					
%% generator cost data																					
%	1	startup	shutdown	n	x1	y1	...	xn	yn												
%	2	startup	shutdown	n	c(n-1)	...	c0														

mpc.gencost = [																					
	2	0	0	2	5	0;																																									
];
																	
																					
																					
																					
																					
																					
																					
																					
																					
