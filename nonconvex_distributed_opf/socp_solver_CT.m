function [cvx_status,sol,exact_max] = socp_solver_CT(net,t0)
%% parameter 
Nbus = net.Nbus;
Nline = net.Nline;
Tpred = net.Tpred;

%% solve socp 
cvx_begin quiet
    %% solver
    % cvx_precision best
    cvx_solver mosek

    %% variable 
    variable s(Nbus,Tpred) complex                % si            
    variable s_gen(Nbus,Tpred) complex 
    variable v(Nbus,Tpred) nonnegative            % vi
    variable S(Nline,Tpred) complex               % Sij
    variable L(Nline,Tpred) nonnegative           % Lij
    
    variable p_bat(Nbus,Tpred) 
    variable e_bat(Nbus,Tpred)
    
    expression S_sum(Nbus,Tpred)                  % sum(Shi - Zhi*Lhi)    
    S_sum(:,:) = 0;    
    for k = 1:Nline
        i = net.fbus(k,1);
        S_sum(i,:) = S_sum(i,:) + S(k,:) - complex(net.r(k),net.x(k))*L(k,:);
    end  
    
    %% objective 
    if net.ncost == 3
        % minimize sum(net.c_gen(:,1)' * (real(s_gen) * net.baseMVA).^2 ) % + sum(net.c_gen(:,2)' * (real(s_gen) * net.baseMVA) )% + sum(net.c_gen(:,2)' * ones(size(real(s_gen)))
        minimize sum(net.c_gen(:,1)' * (real(s_gen)).^2 )
    elseif net.ncost == 2 
        minimize sum(net.c_gen(:,1)' * (real(s_gen) * net.baseMVA)) % + sum(net.c_gen(:,2)' * ones(size(real(s_gen)))
    end
    
    %% constraints
    subject to    
    % power flow constraints
    for k = 1:Nline
        j = net.fbus(k);
        i = net.tbus(k);       
        S(k,:) == s(i,:) + S_sum(i,:);        
        v(i,:) - v(j,:) == 2*real(complex(net.r(k),-net.x(k))*S(k,:)) - (net.r(k)^2 + net.x(k)^2)*L(k,:);
        for tt = 1:Tpred
            norm([2*S(k,tt);L(k,tt)-v(i,tt)],2) <= L(k,tt) + v(i,tt);
        end
    end
    
    0 == s(1,:) + S_sum(1,:);       % s(1) = net power at Bus 1

    % power injection constraints 
    real(s) == real(s_gen) - net.p_demand + p_bat;
    imag(s) == imag(s_gen) - net.q_demand;
    
    repmat(net.p_gen_min,1,Tpred) <= real(s_gen) <= repmat(net.p_gen_max,1,Tpred);
    repmat(net.q_gen_min,1,Tpred) <= imag(s_gen) <= repmat(net.q_gen_max,1,Tpred);
    
    repmat(net.p_bat_min,1,Tpred) <= p_bat <= repmat(net.p_bat_max,1,Tpred);
    repmat(net.e_bat_min,1,Tpred) <= e_bat <= repmat(net.e_bat_max,1,Tpred); 
    e_bat(:,Tpred) == net.e_bat_now;
    if Tpred > 1
        if Tpred == 24
        e_bat(:,24-t0) == net.e_bat_midnight;
        end
        % battery dynamics
        e_bat(:,1) == net.e_bat_now - p_bat(:,1);
        for tt = 2:Tpred
            e_bat(:,tt) == e_bat(:,tt-1) - p_bat(:,tt);
        end
    else
        % error('battery setting is not correct');
    end
    % voltage constraint    
    repmat(net.v_min.^2,1,Tpred) <= v <= repmat(net.v_max.^2,1,Tpred); 
    
cvx_end

sol.obj = cvx_optval;
sol.s = s;
sol.p_gen = real(s_gen);
sol.q_gen = imag(s_gen);
sol.v = v;
sol.p_bat = p_bat;
sol.e_bat = e_bat;
sol.P = zeros(Nbus,Tpred);
sol.Q = zeros(Nbus,Tpred);
sol.L = zeros(Nbus,Tpred);
for iline = 1:net.Nline
    idx = net.tbus(iline);
    sol.L(idx,:) = L(iline,:);
    sol.P(idx,:) = real(S(iline,:));
    sol.Q(idx,:) = imag(S(iline,:));
end
  
%% exactness check 
exact = [];
for k = 1:Nline
    j = net.fbus(k);
    i = net.tbus(k); 
    exact(k,:) = L(k,:).*v(i,:)./abs(S(k,:)).^2;
end    

exact_max = max(exact);

end