function [cvx_status,sol,exact_max] = socp_solver_GT(net,mpc)
%% parameter 
Nbus = net.Nbus;
Nline = net.Nline;

%% solve socp 
cvx_begin quiet
    %% solver 
    cvx_precision best
    cvx_solver mosek
    
    %% variable 
    variable p(Nbus)
    variable q(Nbus)
    variable p_gen(Nbus)             
    variable q_gen(Nbus)
    variable C(Nbus,Nbus) symmetric
    variable S(Nbus,Nbus) skew_symmetric
    
    %% objective 
    if net.ncost == 3
        minimize sum(sum(net.c_gen .* [ (p_gen * mpc.baseMVA).^2 (p_gen * mpc.baseMVA) ones(size(p_gen))]))
    elseif net.ncost == 2
        minimize sum(sum(net.c_gen .* [(p_gen * mpc.baseMVA) ones(size(p_gen))]))
    end
    
    %% constraints
    subject to
    % power flow constraints
    for ii = 1:Nbus
        p(ii) == sum( net.G(ii,:) .* C(ii,:) - net.B(ii,:) .* S(ii,:) );
        q(ii) == sum( - net.B(ii,:) .* C(ii,:) - net.G(ii,:) .* S(ii,:) );
    end
    for ii = 1:Nline
        idxFrom = mpc.branch(ii,1);
        idxTo = mpc.branch(ii,2); 
        norm([2*C(idxFrom,idxTo); 2*S(idxFrom,idxTo); C(idxFrom,idxFrom)-C(idxTo,idxTo) ],2) <= C(idxFrom,idxFrom)+C(idxTo,idxTo);
    end
    
    % power injection constraints     
    p == p_gen - net.p_demand;
    q == q_gen - net.q_demand;    
    net.p_gen_min <= p_gen <= net.p_gen_max;
    net.q_gen_min <= q_gen <= net.q_gen_max;

    % voltage constraint
    net.v_min.^2 <= diag(C) <= net.v_max.^2;
    
cvx_end

sol.obj = cvx_optval;
sol.p = p;
sol.q = q;
sol.p_gen = p_gen;           
sol.q_gen = q_gen;
sol.C = C;
sol.S = S;


%% exactness check 
exact = [];
for ii = 1:Nline
    idxFrom = mpc.branch(ii,1);
    idxTo = mpc.branch(ii,2); 
    exact(ii,1) = C(idxFrom,idxFrom)*C(idxTo,idxTo)/(C(idxFrom,idxTo)^2+S(idxFrom,idxTo)^2);
end    

exact_max = max(exact);

end