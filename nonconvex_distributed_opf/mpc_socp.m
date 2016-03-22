function sol_socp = mpc_socp(net,p_pred,q_pred)

% closed loop operation with convex relaxation approach 

net.e_bat_now = net.e_bat_midnight;

for ii = 1:net.Ncalc 

    net.p_demand = p_pred{ii};
    net.q_demand = q_pred{ii};
    t0 = ii-1;
    
    [cvx_status,sol_CT,exact_max] = socp_solver_CT(net,t0);
    disp('- SOCP (CalT formulation):')
    fprintf('Objective Value = %.2f\n',sol_CT.obj);
    fprintf('Exactness Worstcase = %.6f\n',mean(exact_max));
    
    sol_socp{ii} = sol_CT;
    
    % parameter update
    net.e_bat_now = sol_CT.e_bat(:,1);

end

end