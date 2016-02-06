close all;
clear all;
clc;

%% load network data from GT
% mpc = case2Tree;
mpc = case3Tree;

%% pre-process network

net = pre_opf_net(mpc);

lambda = 0.9; %0.81;
net.p_demand = net.p_demand * lambda;
net.q_demand = net.q_demand * lambda;
fprintf('lambda = %.2f\n',lambda);
% return;
%% OPF-SOCP from CalT 

[cvx_status,sol_CT,exact_max] = socp_solver_CT(net);
disp('- SOCP (CalT formulation):')
fprintf('Objective Value = %.2f\n',sol_CT.obj);
fprintf('Exactness Worstcase = %.4f\n',exact_max);

%% OPF non-convex distributed BCD algorithm 
% return;
%distr_admm_solver(net)
return;
[sol_bcd,dist_x,dist_z,dist_dual] = distr_opf_solver(net);

