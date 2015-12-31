function [cvx_status,sol,exact_max] = socp_solver_CT(net,mpc)
%% parameter 
Nbus = net.Nbus;
Nline = net.Nline;

%% solve socp 
cvx_begin quiet
    %% solver
    cvx_precision best
    cvx_solver mosek

    %% variable 
    variable s(Nbus) complex                % si            
    variable s_gen(Nbus) complex 
    variable v(Nbus) nonnegative            % vi
    variable S(Nline) complex               % Sij
    variable L(Nline) nonnegative           % Lij
    
    expression S_sum(Nbus)                  % sum(Shi - Zhi*Lhi)    
    S_sum(:) = 0;    
    for k = 1:Nline
        i = mpc.branch(k,1);
        S_sum(i) = S_sum(i) + S(k) - complex(net.r(k),net.x(k))*L(k);
    end  
    
    %% objective 
    if net.ncost == 3
        minimize sum(sum(net.c_gen .* [(real(s_gen) * mpc.baseMVA).^2 (real(s_gen) * mpc.baseMVA) ones(size(real(s_gen)))]))
    elseif net.ncost == 2
        minimize sum(sum(net.c_gen .* [(real(s_gen) * mpc.baseMVA) ones(size(real(s_gen)))]))
    end
    
    %% constraints
    subject to    
    % power flow constraints
    for k = 1:Nline
        j = mpc.branch(k,1);
        i = mpc.branch(k,2);       
        S(k) == s(i) + S_sum(i);        
        v(i) - v(j) == 2*real(complex(net.r(k),-net.x(k))*S(k)) - (net.r(k)^2 + net.x(k)^2)*L(k);
        norm([2*S(k);L(k)-v(i)],2) <= L(k) + v(i);
    end
    
    0 == s(1) + S_sum(1);       % s(1) = net power at Bus 1

    % power injection constraints 
    s == s_gen - complex(net.p_demand,net.q_demand);

    net.p_gen_min <= real(s_gen) <= net.p_gen_max;
    net.q_gen_min <= imag(s_gen) <= net.q_gen_max;
    
    % voltage constraint    
    net.v_min.^2 <= v <= net.v_max.^2; 
    
cvx_end

sol.obj = cvx_optval;
sol.s = s;
sol.s_gen = s_gen;
sol.v = v;
sol.S = S;
sol.L = L;

%% exactness check 
exact = [];
for k = 1:Nline
    j = mpc.branch(k,1);
    i = mpc.branch(k,2); 
    exact(k,1) = L(k)*v(i)/abs(S(k))^2;
end    

exact_max = max(exact);

end