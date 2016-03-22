% :: main script for real-time distributed opf 

close all;
clear all;
clc;

%% load network data from GT
% mpc = case2Tree;
% mpc = case3Tree;
mpc = case9Tree;

%% pre-process network

net = pre_opf_net(mpc);

net.Tpred = 24;
fprintf('Prediction step = %d \n',net.Tpred);
p_factor = stat_hr_profile(net.Nbus,net.Tpred);
q_factor = ones(net.Nbus,net.Tpred); 

net.p_demand = repmat(net.p_demand,1,net.Tpred) .* p_factor;
net.q_demand = repmat(net.q_demand,1,net.Tpred) .* q_factor;
fprintf('Network prepared.\n');

%% demand prediction 

net.Ncalc = 24;
PERTMAX = 0.04;
rng('default');     % generate same perturbation seed 
[p_pred,q_pred] = loadpred(net,PERTMAX);

%% socp solver 
sol_socp = mpc_socp(net,p_pred,q_pred);

net.p_gen_ext = zeros(net.Nbus,net.Tpred);
for ii = 1:net.Tpred
    net.p_gen_ext(1,ii) = sol_socp{1}.p_gen(1,ii);
end

%% opf solver 
sol_opf_real_time = distr_opf_solver(net,sol_socp,p_pred,q_pred);

%% offline plot
ii = 1;
iterplot(sol_opf_real_time{ii},sol_socp{ii});

%% closed loop plot

net_mpc.e_bat_now = net.e_bat_midnight;
net_mpc.ri = net.ri;
net_mpc.Tpred = net.Ncalc;
net_mpc.baseMVA = net.baseMVA;
t0 = 0;

for ii = 1:net.Ncalc
    sol_socp_mpc.e_bat(:,ii) = sol_socp{ii}.e_bat(:,1);
    sol_socp_mpc.p_gen(:,ii) = sol_socp{ii}.p_gen(:,1);
    sol_socp_mpc.L(:,ii) = sol_socp{ii}.L(:,1);
    sol_socp_mpc.p_bat(:,ii) = sol_socp{ii}.p_bat(:,1);
    sol_socp_mpc.obj(1,ii) = sol_socp{ii}.obj(end);
    sol_socp_mpc.res(1,ii) = norm( sol_socp{ii}.v .* sol_socp{ii}.L - sol_socp{ii}.P.^2 - sol_socp{ii}.Q.^2 );
    net_mpc.p_demand(:,ii) = p_pred{ii}(:,1);
end    

sol_disp = sol_opf_real_time;
for ii = 1:net.Ncalc
    sol_disp_mpc.e_bat(:,ii) = sol_disp{ii}.e_bat(:,1);
    sol_disp_mpc.p_gen(:,ii) = sol_disp{ii}.p_gen(:,1);
    sol_disp_mpc.L(:,ii) = sol_disp{ii}.L(:,1);
    sol_disp_mpc.p_bat(:,ii) = sol_disp{ii}.p_bat(:,1);
    if ii > 1
        sol_disp_mpc.res(1,ii) = sqrt(sol_disp{ii}.resPrim(end)^2 + sol_disp{ii}.resPF(end)^2);
    else
        sol_disp_mpc.res(1,ii) = 1e-4;
    end
    sol_disp_mpc.obj(1,ii) = sol_disp{ii}.obj(end);
    net_mpc.p_demand(:,ii) = p_pred{ii}(:,1);
end
% simplot(net_mpc,sol_disp_mpc,t0)
% sol_socp_mpc.e_bat - sol_opf_mpc.e_bat

offline_mpc.p_demand = p_pred{1};
offline_mpc.e_bat = sol_disp{1}.e_bat;
offline_mpc.p_gen = sol_disp{1}.p_gen;
rtplot(offline_mpc,net_mpc,sol_disp_mpc,sol_socp_mpc)
