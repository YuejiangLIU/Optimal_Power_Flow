function ctrl = centralized_controller(mpc,b0,p_load_rated,q_load_rated,t0,bnow,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dev_V = mpc.dev_V;
dev_Vmin = dev_V(1);
dev_Vmax = dev_V(2);
p_generator_lowbound = mpc.p_generator_lowbound;
s_generator_upbound = mpc.s_generator_upbound;
p_capacitor_rated = mpc.p_capacitor_rated;
q_capacitor_upbound = mpc.q_capacitor_upbound;
if mode > 1
    B_upbound = mpc.B_upbound;
end

% >>>> Version Adaption <<<< 
R_Line = mpc.branch(:,3)/mpc.base_Z;
X_Line = mpc.branch(:,4)/mpc.base_Z;

yalmip('clear');
IterLim = Inf;
clear S_sum;
[Nbus,T] = size(s_generator_upbound);
Nline = Nbus -1;

% *********************** yalmip ************************************ 
% --- variables ---
s = sdpvar(Nbus,T,'full','complex');
s_load = sdpvar(Nbus,T,'full','complex');
s_capacitor = sdpvar(Nbus,T,'full','complex');
s_generator = sdpvar(Nbus,T,'full','complex');
if mode > 1
    s_storage = sdpvar(Nbus,T,'full');
    b = sdpvar(Nbus,T+1,'full');
end
v = sdpvar(Nbus,T,'full');
S = sdpvar(Nline,T,'full','complex');
L = sdpvar(Nline,T,'full');

% --- expressions ---  
for k = 1:Nline
    i = mpc.branch(k,1);
    h = mpc.branch(k,2);
    if ~exist('S_sum','var') || size(S_sum,1) < i
        S_sum(i,:) = S(k,:) - complex(R_Line(k),X_Line(k)).*L(k,:);
    else
        S_sum(i,:) = S_sum(i,:) + S(k,:) - complex(R_Line(k),X_Line(k)).*L(k,:);
    end
end  

% --- objective ---
objective = 0;
c1 = 1; 
c2 = 2;
% objective = c1 * sum(sum( real(s_generator).*real(s_generator) )) + c2 * sum( real(s(1,:)).*real(s(1,:)) );
% objective = objective + 1e5 * sum(sum(real(s)));
% objective = objective + c2 * sum(real(s(1,:)));
objective = objective + sum(sum(real(s)));

% --- constraints ---
constraints = [];

% constraints = constraints + [real(s(1,:))>=0];

% power flow constraints
for tt = 1:T
    for k = 1:Nline
        j = mpc.branch(k,1);
        i = mpc.branch(k,2);
        if i > size(S_sum,1)
            constraints = constraints + [S(k,:) == s(i,:)];
        else
            constraints = constraints + [S(k,:) == s(i,:) + S_sum(i,:)];
        end
        constraints = constraints + [v(i,:) - v(j,:) == 2*real(complex(R_Line(k),-X_Line(k))*S(k,:)) - (R_Line(k)^2 + X_Line(k)^2)*L(k,:)];
        constraints = constraints + [norm([2*S(k,tt);L(k,tt)-v(i,tt)],2) <= L(k,tt) + v(i,tt)];
    end
end
constraints = constraints + [ s(1,:) + S_sum(1,:) == 0];       % s(1) = net power at Bus 1
% power injection constraints
if mode > 1
    constraints = constraints + [s(2:end,:) == s_load(2:end,:) + s_capacitor(2:end,:) + s_generator(2:end,:) + s_storage(2:end,:)];
elseif mode == 1
    constraints = constraints + [s(2:end,:) == s_load(2:end,:) + s_capacitor(2:end,:) + s_generator(2:end,:)];
end

constraints = constraints + [real(s_load) == p_load_rated, imag(s_load) == q_load_rated];
constraints = constraints + [real(s_capacitor) == p_capacitor_rated];
for i = 1:Nbus
    if q_capacitor_upbound(i,:) > 0 
        constraints = constraints + [0 <= imag(s_capacitor(i,:)) <= q_capacitor_upbound(i,:)];
    else 
        constraints = constraints + [imag(s_capacitor(i,:)) == 0];
    end
    if s_generator_upbound(i,:) > 0 
        constraints = constraints + [real(s_generator(i,:)) >= p_generator_lowbound(i,:)];
        constraints = constraints + [abs(s_generator(i,:)) <= s_generator_upbound(i,:)];
    else 
        constraints = constraints + [s_generator(i,:) == 0];
    end
    if mode > 1
    if B_upbound(i,:) > 0
        constraints = constraints + [- B_upbound(i,:) / 4 <= s_storage(i,:) <= B_upbound(i,:) / 4];
        constraints = constraints + [0 <= b(i,2:end) <= B_upbound(i,:)];
    else
        constraints = constraints + [s_storage(i,:) == 0];
    end
    end
end
% battery constraints
if mode > 1
constraints = constraints + [b(:,1) == bnow];
%constraints = constraints + [b(:,T+1-t0) == b0];
%constraints = constraints + [b(:,1) == b(:,T+1)];
for tt = 1:T
    constraints = constraints + [b(:,tt+1) == b(:,tt) - s_storage(:,tt)];
end
end
% voltage constraint    
constraints = constraints + [v(1,:) == 1];
constraints = constraints + [((1-dev_Vmin))^2 <= v(2:end,:) <= ((1+dev_Vmax))^2]; 

% --- options ---
%ops = sdpsettings('solver','sedumi','sedumi.eps',1e-8);
%ops = sdpsettings('solver','sdpt3','verbose',1);
ops = sdpsettings('solver','sedumi','verbose',1);
%ops = sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.OptimalityTol',1e-9,'gurobi.BarIterLimit',IterLim,'verbose',1);
%ops = sdpsettings('solver','gurobi','verbose',1);

% --- solve --- 
diagnostics = optimize(constraints, objective, ops);

if diagnostics.problem == 0
 %disp('Feasible')
elseif diagnostics.problem == 1
 disp('Solver met problem. Please check objective & constrains. ');
 error('Infeasible')
else
 diagnostics.problem
 disp('Solver met problem. Please check objective & constrains. ');    
 error('Something else happened')
end

%ctrl.b = value(b);
ctrl.obj = value(objective);
ctrl.L = value(L);
ctrl.v = value(v);
ctrl.s = value(s);
ctrl.gen = sum(value(s_generator))+value(s(1,:));
ctrl.S = value(S);

ctrl.exact = ctrl.L .* ctrl.v(2:end,:) ./ abs(ctrl.S).^2;

double(objective)
loss = sum(double(sum(L.*repmat(R_Line,1,T))))        % verification for objective function


end