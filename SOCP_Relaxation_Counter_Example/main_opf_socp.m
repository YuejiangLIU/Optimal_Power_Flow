close all;
clear all;
clc;

%% load network data from GT
mpc = case2Tree;
% mpc = case3Tree;

%% pre-process network

net = pre_opf_net(mpc);

lambda = 0.81;
net.p_demand = net.p_demand * lambda;
net.q_demand = net.q_demand * lambda;
fprintf('lambda = %.2f\n',lambda);

%% OPF-SOCP from GT 
disp('============= computation =============')

[cvx_status,sol_GT,exact_max] = socp_solver_GT(net,mpc);
disp('- SOCP (GT formulation):')
fprintf('Objective Value = %.2f\n',sol_GT.obj);
fprintf('Exactness Worstcase = %.4f\n',exact_max);

%% OPF-SOCP from CalT 

[cvx_status,sol_CT,exact_max] = socp_solver_CT(net,mpc);
disp('- SOCP (CalT formulation):')
fprintf('Objective Value = %.2f\n',sol_CT.obj);
fprintf('Exactness Worstcase = %.4f\n',exact_max);

% pause;
%% comparison
format shortG
disp('============= comparison =============')
disp('form:    CalT           GT');
disp('voltage:');
disp([sqrt(sol_CT.v) sqrt(diag(sol_GT.C))]);
disp('p_gen:');
disp([real(sol_CT.s_gen) sol_GT.p_gen]);
disp('q_gen:');
disp([imag(sol_CT.s_gen) sol_GT.q_gen]);
disp('p:');
disp([real(sol_CT.s) sol_GT.p]);
