function mpc = case9Tree	
%BASED ON																				
%CASE9    Power flow data for 9 bus, 3 generator case.																					
%   Please see CASEFORMAT for details on the case file format.																					
%																					
%   Based on data from Joe H. Chow's book, p. 70.																					
																					
%   MATPOWER																					
%   $Id: case9.m,v 1.11 2010/03/10 18:08:14 ray Exp $																					
																					
%% MATPOWER Case Format : Version 2																					
mpc.version = '2';																					
																					
%%-----  Power Flow Data  -----%%																					
%% system MVA base																					
mpc.baseMVA = 100;																					
																					
%% bus data																					
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin								
mpc.bus = [																					
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	2	2	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	3	2	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	4	1	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	5	1	90	30	0	0	1	1	0	345	1	1.1	0.9;								
	6	1	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	7	1	100	35	0	0	1	1	0	345	1	1.1	0.9;								
	8	1	0	0	0	0	1	1	0	345	1	1.1	0.9;								
	9	1	125	50	0	0	1	1	0	345	1	1.1	0.9;								
];																					
																					
%% generator data																					
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [																					
	1	0	0	300	100	1	100	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	2	163	0	300	100	1	100	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
	3	85	0	300	100	1	100	1	270	10	0	0	0	0	0	0	0	0	0	0	0;
];																					
																					
%% branch data																					
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax								
mpc.branch = [																					
	1	4	0	0.0576	0	9900	250	250	0	0	1	-360	360;								
	4	5	0.017	0.092	0.158	9900	250	250	0	0	1	-360	360;								
	5	6	0.039	0.17	0.358	9900	150	150	0	0	1	-360	360;								
	6	3	0	0.0586	0	9900	300	300	0	0	1	-360	360;								
	8	7	0.0085	0.072	0.149	9900	250	250	0	0	1	-360	360;								
	8	2	0	0.0625	0	9900	250	250	0	0	1	-360	360;								
	9	8	0.032	0.161	0.306	9900	250	250	0	0	1	-360	360;								
	4	9	0.01	0.085	0.176	9900	250	250	0	0	1	-360	360;								
];																					
																					
%%-----  OPF Data  -----%%																					
%% area data																					
%	area	refbus																			
mpc.areas = [																					
	1	5;																			
];																					
																					
%% generator cost data																					
%	1	startup	shutdown	n	x1	y1	...	xn	yn												
%	2	startup	shutdown	n	c(n-1)	...	c0														
%{
mpc.gencost = [																					
	2	1500	0	2	5	0;														
	2	2000	0	2	1.2	0;														
	2	3000	0	2	1	0;														
];																					
%}
mpc.gencost = [																					
	2	1500	0	3	0.5     0   0;														
	2	2000	0	3	0.12	0   0;														
	2	3000	0	3	0.1     0   0;	
];

%% battery data																																																														
%	bus	Pmax	Pmin	Emax    Emin    E0      
mpc.bat = [																					
	3    20      -20     60      0     20           
];	    
