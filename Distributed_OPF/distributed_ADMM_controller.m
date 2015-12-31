function sol = distributed_ADMM_controller(load,B0,t0,bnow,mpc,rho,Niter,objTarget)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ADMM-based distributed opf controller 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yalmip('clear');

[Nbus,T] = size(load);

% Rated Load
p_load_rated = - load*mpc.pf;
q_load_rated = - load*sqrt(1-mpc.pf^2);

% Initialization 
x = mpc.x;
z = mpc.z;
lamda = mpc.lamda;

%waitname = strcat('Prediction at T = ',num2str(t0));
%h = waitbar(0,waitname);
Ncurrentiter = 0;
%figure;
tStart = tic; 
format shortG;


for niter = Ncurrentiter + 1:Ncurrentiter + Niter 
%% ADMM x-update 

for ii = 1:Nbus
    
    Na = mpc.A(ii,1);
    Nc = mpc.C(ii,1);
    
    if niter == 1
% vvvvvvvvvvvvvvvvvvvvvvvv YALMIP Optimizer Start vvvvvvvvvvvvvvvvvvvvvvv    
    % ---- parameters ---- 
    clear lamdaAA lamdaBB lamdaCC ZZ;
    lamdaBB = sdpvar(4,1,T,'full');
    if Na > 0 
        lamdaAA = sdpvar(1,Na,T,'full');
    else 
        clear lamdaAA;
    end
    if Nc > 0 
        lamdaCC = sdpvar(3,Nc,T,'full');
    else 
        clear lamdaCC;
    end
    ZZ = sdpvar(4,Nbus,T,'full');
    
    % ---- variable ----           
    clear xB;
    xB = sdpvar(6,1,T,'full');    % [v; L; P; Q; p; q]
    
    clear xA;
    if Na>0
        xA = sdpvar(1,1,T,'full');             % [va]
    else
        xA = [];
    end
    
    clear xC;
    if Nc>0
        xC = sdpvar(3,Nc,T,'full');       % [Lc; Pc; Qc]
    else 
        xC = [];
    end
    
    BAT = sdpvar(1,T+1,'full');
    p_battery = sdpvar(1,T,'full');
    
    % ---- auxiliary variable ----
    p_generator = sdpvar(1,T,'full');
    q_generator = sdpvar(1,T,'full');
    q_capacitor = sdpvar(1,T,'full');
    
    
    % ---- expression in function and constraint ----
    % >>> sub objective <<<
    clear fi;
    %fi = norm(p_generator,2);               % flatten generation variance 
    fi = sum_square(p_generator); 
    
    clear productB;
    productB = sum(sum(lamdaBB.*xB(1:4,1,:))); 
    
    clear productA;
    if Na > 0     
        productA = sum(lamdaAA.*xA);
    else
        productA = 0;
    end
    
    clear productC;
    if Nc > 0
        if Nc == 1
            productC = sum(sum(lamdaCC.*xC));
        else        % Nc > 1 => xC is a matrix 
            productC = sum(sum(sum(lamdaCC.*xC)));
        end
    else
        productC = 0;
    end
    
    clear sumsqB;
    sumsqB = sum(sum_square(xB(1:4,1,:) - ZZ(:,ii,:)));
    
    clear sumsqA;
    if Na > 0 
        indA = mpc.A(ii,2);
        sumsqA = sum_square( xA - ZZ(1,indA,:) );
    else
        sumsqA = 0;
    end
    
    clear sumsqC;
    sumsqC = 0;
    if Nc > 0 
        for cc = 1:Nc               % enumerate chilren 
            indC = mpc.C(ii,cc+1);      % who is the child
            sumsqC = sumsqC + sum(sum_square( xC(:,cc,:) - ZZ(2:4,indC,:)));
        end
    else 
        sumsqC = 0;
    end
            
    clear sumPC sumQC;
    sumPC = 0;
    sumQC = 0;
    if Nc > 0 
        for cc = 1:Nc
            indC = mpc.C(ii,cc+1);
            sumPC = sumPC + ( xC(2,cc,:) - xC(1,cc,:) * mpc.R(indC) );
            sumQC = sumQC + ( xC(3,cc,:) - xC(1,cc,:) * mpc.X(indC) );
        end
    else 
        sumPC = 0;
        sumQC = 0;
    end
           
    % ---- objective ----
    objective = fi + productB + productA + productC ...
        + rho/2 * ( sumsqB + sumsqA + sumsqC ); 
    
    % ---- constraint ----
    constraints = [];
    % power flow constraints
    if Na > 0
        constraints = constraints + [xA == xB(1,1,:) - 2 * ( mpc.R(ii) * xB(3,1,:) + mpc.X(ii) * xB(4,1,:) ) + xB(2,1,:) * ( mpc.R(ii)^2 + mpc.X(ii)^2 )];
    else 
        constraints = constraints + [xB(1,1,:) == 1];
        constraints = constraints + [xB(2,1,:) == 0];
        constraints = constraints + [xB(3,1,:) == 0];
        constraints = constraints + [xB(4,1,:) == 0];
    end        
    constraints = constraints + [sumPC + xB(5,1,:) == xB(3,1,:)];
    constraints = constraints + [sumQC + xB(6,1,:) == xB(4,1,:)];
    for tt = 1:T
    constraints = constraints + [norm([ 2*xB(3,1,tt); 2*xB(4,1,tt); xB(2,1,tt)-xB(1,1,tt) ],2) <= xB(2,1,tt)+xB(1,1,tt)];
    end
    % voltage constraints    
    if ii > 1 
        constraints = constraints + [(1-mpc.dev_V)^2 <= xB(1,1,:) <= (1+mpc.dev_V)^2];
    else 
        constraints = constraints + [xB(1,1,:) == 1];
    end    
    if Na > 0 
        indA = mpc.A(ii,2);
        if indA > 1
            constraints = constraints + [(1-mpc.dev_V)^2 <= xA(1,1,:) <= (1+mpc.dev_V)^2];
        else
            constraints = constraints + [xA(1,1,:) == 1];
        end
    end
    % power injection constraints     
    for tt = 1:T 
        constraints = constraints + [xB(5,1,tt) == p_load_rated(ii,tt) + p_generator(1,tt) + p_battery(1,tt)];
        constraints = constraints + [xB(6,1,tt) == q_load_rated(ii,tt) + q_generator(1,tt) + q_capacitor(1,tt)];
    end    
    if mpc.q_capacitor_upbound(ii) > 0 
        constraints = constraints + [0 <= q_capacitor <= mpc.q_capacitor_upbound(ii,:)];
    else 
        constraints = constraints + [q_capacitor == 0];
    end
    if mpc.s_generator_upbound(ii) > 0 
        constraints = constraints + [p_generator >= mpc.p_generator_lowbound(ii,:)];
        for tt = 1:T
            constraints = constraints + [norm([p_generator(1,tt);q_generator(1,tt)],2) <= mpc.s_generator_upbound(ii,tt)];
        end
    else 
        constraints = constraints + [p_generator == 0];
        constraints = constraints + [q_generator == 0];
    end
    % battery constraints
    constraints = constraints + [- mpc.B_upbound / 4 <= p_battery <= mpc.B_upbound / 4];
    constraints = constraints + [0 <= BAT <= mpc.B_upbound];
    constraints = constraints + [BAT(1,T+1-t0) == B0(ii)];
    constraints = constraints + [BAT(1,1) == bnow(ii)];
    constraints = constraints + [BAT(1,1) == BAT(1,T+1)];
    for tt = 1:T
        constraints = constraints + [BAT(1,tt+1) == BAT(1,tt) - p_battery(1,tt)];
    end
    
    % ---- settings ----
    ops = sdpsettings('solver','gurobi');
    %ops = sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.OptimalityTol',1e-9);%,'verbose',0);
    %ops = sdpsettings('solver','mosek');%,'verbose',0);

    % ---- modeling ----
    if Na < 1
        xoptimizer{ii} = optimizer(constraints, objective, ops,...
            {lamdaBB,lamdaCC,ZZ},{xB,xC,fi,BAT,p_generator,q_generator});        
    elseif Nc < 1
        xoptimizer{ii} = optimizer(constraints, objective, ops,...
            {lamdaAA,lamdaBB,ZZ},{xA,xB,fi,BAT,p_generator,q_generator});        
    else
        xoptimizer{ii} = optimizer(constraints, objective, ops,...
            {lamdaAA,lamdaBB,lamdaCC,ZZ},{xA,xB,xC,fi,BAT,p_generator,q_generator});
    end
    
    
% ^^^^^^^^^^^^^^^^^^^^^^ YALMIP Optimizer End ^^^^^^^^^^^^^^^^^^^^^^^^^
    end

    % --- YALMIP Controller --- 
    if Na < 1
        xctrl = xoptimizer{ii};
        [xopt, errcode] = xctrl{{lamda(ii).B,lamda(ii).C,z}};
        xB = xopt{1};
        xC = xopt{2};
        fi = xopt{3};
        battery = xopt{4};
        p_generator = xopt{5};
        q_generator = xopt{6};
    elseif Nc <1 
        xctrl = xoptimizer{ii};
        [xopt, errcode] = xctrl{{lamda(ii).A,lamda(ii).B,z}};
        xA = xopt{1};
        xB = xopt{2};
        fi = xopt{3};
        battery = xopt{4};
        p_generator = xopt{5};
        q_generator = xopt{6};
    else 
        xctrl = xoptimizer{ii};
        [xopt, errcode] = xctrl{{lamda(ii).A,lamda(ii).B,lamda(ii).C,z}};
        xA = xopt{1};
        xB = xopt{2};
        xC = xopt{3};
        fi = xopt{4};
        battery = xopt{5};
        p_generator = xopt{6};
        q_generator = xopt{7};
    end 
    clear xctrl;

    if errcode ~= 0
        disp('***** Gurobi warning. *****');
        disp(yalmiperror(errcode));
        fprintf('niter = %d     ii = %d\n',niter,ii);
        disp('lamda(ii).A =');
        lamda(ii).A
        disp('lamda(ii).B =');
        lamda(ii).B
        disp('lamda(ii).C =');
        lamda(ii).C
        disp('***************************');
    end
    
x(ii).B = double(xB);
if Na>0 
    x(ii).A = double(xA);
end
if Nc>0 
    x(ii).C = double(xC);
end

objLocal(ii,1) = double(fi);
batLocal(ii,:) = double(battery);
p_gen(ii,:) = p_generator;
q_gen(ii,:) = q_generator;

% --- result check ---
%{
z
disp('lamda.B');
lamda.B
disp('x.B');
x.B
if Na>0 
    disp('lamda.A');
    lamda.A 
    disp('x.A');
    x.A 
end
if Nc > 0 
    disp('lamda.C');
    lamda.C
    disp('x.C');
    x.C
end
%}

end;

objTotal(niter,1) = sum(objLocal);
batTotal = sum(batLocal);
%% ADMM z-update 
zold = z;
for ii = 1:Nbus
    
    Na = mpc.A(ii,1);
    Nc = mpc.C(ii,1);

    % --- update z_v ---
    if Nc > 0
        facA = rho/2 * (1 + Nc);
        facB = - rho*x(ii).B(1,1,:) - lamda(ii).B(1,1,:);
        for cc = 1:Nc
            indC = mpc.C(ii,cc+1);
            facB = facB - rho*x(indC).A(1,1,:) - lamda(indC).A(1,1,:);
        end
        z(1,ii,:) = -facB/facA/2;
    else 
        z(1,ii,:) = x(ii).B(1,1,:) + lamda(ii).B(1,1,:)/rho;
    end

    clear facA;
    clear facB;

    % --- update z_L, z_P, z_Q ---
    if Na > 0
        indA = mpc.A(ii,2);         % who is ii's ancestor 
        nbC = find(mpc.C(indA,:) == ii,1,'last') - 1; % ii is his ancestor's No. # child 
        facA = rho;
        facB(1,1,:) = rho*(-x(ii).B(2,1,:) - x(indA).C(1,nbC,:))-lamda(ii).B(2,1,:) - lamda(indA).C(1,nbC,:);
        facB(2,1,:) = rho*(-x(ii).B(3,1,:) - x(indA).C(2,nbC,:))-lamda(ii).B(3,1,:) - lamda(indA).C(2,nbC,:);
        facB(3,1,:) = rho*(-x(ii).B(4,1,:) - x(indA).C(3,nbC,:))-lamda(ii).B(4,1,:) - lamda(indA).C(3,nbC,:);
        z(2,ii,:) = -facB(1,1,:)/facA/2;
        z(3,ii,:) = -facB(2,1,:)/facA/2;
        z(4,ii,:) = -facB(3,1,:)/facA/2;
    else
        z(2:4,ii,:) = x(ii).B(2:4,1,:) + lamda(ii).B(2:4,1,:)/rho;
    end
    clear facA;
    clear facB;

% --- result check ---
%{
z
disp('lamda.B');
lamda.B
disp('x.B');
x.B
if Na>0 
    disp('lamda.A');
    lamda.A 
    disp('x.A');
    x.A 
end
if Nc > 0 
    disp('lamda.C');
    lamda.C
    disp('x.C');
    x.C
end
%}
end


%% Multiplier update

for ii = 1:Nbus
    
    Na = mpc.A(ii,1);
    Nc = mpc.C(ii,1);
        
    lamda(ii).B = lamda(ii).B + rho * ( x(ii).B(1:4,1,:) - z(:,ii,:) );
    
    if Na > 0
        indA = mpc.A(ii,2);
        lamda(ii).A = lamda(ii).A + rho * ( x(ii).A - z(1,indA,:) );
    end
    
    if Nc > 0
        for cc = 1:Nc
            indC = mpc.C(ii,cc+1);
            lamda(ii).C(:,cc,:) = lamda(ii).C(:,cc,:) + rho * ( x(ii).C(:,cc,:) - z(2:4,indC,:) );
        end
    end

end

% --- result check ---
%{
z
disp('lamda.B');
lamda.B
disp('x.B');
x.B
if Na>0 
    disp('lamda.A');
    lamda.A 
    disp('x.A');
    x.A 
end
if Nc > 0 
    disp('lamda.C');
    lamda.C
    disp('x.C');
    x.C
end
%}

% --- residual ---
res.B(niter,1) = 0;
for ii = 1:Nbus
    %res.B(niter,1) = res.B(niter,1) + norm( x(ii).B-z(:,ii,:),2 ); 
    res.B(niter,1) = res.B(niter,1) + sum(sum_square(x(ii).B(1:4,1,:)-z(:,ii,:)));
end

res.A(niter,1) = 0;
res.C(niter,1) = 0;
for ii = 1:Nbus
    Na = mpc.A(ii,1);
    if Na > 0
        indA = mpc.A(ii,2);
        res.A(niter,1) = res.A(niter,1) + sum_square( x(ii).A - z(1,indA,:)); 
    end
    
    Nc = mpc.C(ii,1);
    if Nc > 0
        for cc = 1:Nc
            indC = mpc.C(ii,cc+1);
            res.C(niter,1) = res.C(niter,1) + sum(sum_square( x(ii).C(:,cc,:) - z(2:4,indC,:)));
        end
    end
end

resPrimal = sqrt(res.A(niter) + res.B(niter) + res.C(niter));
res.prim(niter,1) = resPrimal;
resDual = rho*sqrt(sum(sum(sum_square(z-zold))));
res.dual(niter,1) = resDual;

% --- iteraction process --- 
%{
xx = 1:niter;
clf;
plot(xx,res.prim,'r-',xx,res.dual,'g-');
set(gca,'yscale','log');
legend('Primal Residual',...
    'Dual Residual','Location','southwest');
grid on;
xlabel('Iteration');
%}

% --- adaptive rho ---
%{
mu = 10;
tau = 2;
if rho < 1e2
    if res.prim(niter) > mu*res.dual(niter)
        rho = rho*tau;
    elseif res.dual(niter) > mu*res.prim(niter)
        rho = rho/tau;
    end
end
%}

%waitbar(niter/Niter+Ncurrentiter);

end

Ncurrentiter = niter;

tADMM = toc(tStart);
loopTime = tADMM/Niter;

%close(h);
%close all;


%% Post Process

for ii = 1:Nbus    
x_v(ii,:) = x(ii).B(1,1,:);
x_L(ii,:) = x(ii).B(2,1,:);
x_P(ii,:) = x(ii).B(3,1,:);
x_Q(ii,:) = x(ii).B(4,1,:);
x_p(ii,:) = x(ii).B(5,1,:);
x_q(ii,:) = x(ii).B(6,1,:);
end

z_v = z(1,:,:);
z_L = z(2,:,:);
z_P = z(3,:,:);
z_Q = z(4,:,:);

%{
a_v = zeros(Nbus,1);
c_L = zeros(Nbus,1);
c_P = zeros(Nbus,1);
c_Q = zeros(Nbus,1);

for ii = 1:Nbus
    Na = A(ii,1);    
    if Na>0
        indA = A(ii,2);
        nbC = find(C(indA,:) == ii,1,'last') - 1; 
        a_v(indA,nbC) = x(ii).A;
    end
    
    Nc = C(ii,1);
    if Nc> 0
        for cc = 1:Nc
            indC = C(ii,cc+1);
            c_L(indC,1) = x(ii).C(1,cc);
            c_P(indC,1) = x(ii).C(2,cc);
            c_Q(indC,1) = x(ii).C(3,cc);
        end
    end
    
end
   
disp('v ->'); 
%[x_v z_v x_v./z_v]
[x_v z_v a_v]
disp('L ->'); 
%[x_L z_L x_L./z_L]
[x_L z_L c_L]
disp('P ->'); 
%[x_P z_P x_P./z_P]
[x_P z_P c_P]
disp('Q ->'); 
%[x_Q z_Q x_Q./z_Q]
[x_Q z_Q c_Q]
disp('p ->'); 
% [x_p z_p x_p./z_p]
[x_p z_p]
disp('q ->'); 
%[x_q z_q x_q./z_q]
[x_q z_q]
%} 
% --- Exact Check ---
P2Q2_x = x_P.^2 + x_Q.^2;
vL_x = x_v.*x_L;
x_exact = vL_x./P2Q2_x;
P2Q2_z = z_P.^2 + z_Q.^2;
vL_z = z_v.*z_L;
z_exact = vL_z./P2Q2_z;

%x_check = table(P2Q2_x,vL_x,x_exact)
%z_check = table(P2Q2_z,vL_z,z_exact)

x_check = [min(min(x_exact)) max(max(x_exact))];
z_check = [min(min(z_exact)) max(max(z_exact))];

%p_station = x_p(1,:) - p_load_rated(1,:);
p_station = p_gen(1,:);
q_station = q_gen(1,:);

%% distributed prediction simulation 
figPred = figure;
t = t0:1:t0+T;
plot(t,[p_station p_station(1)],'--d',...
    t,-[sum(p_load_rated) sum(p_load_rated(:,1))],'--*',...
    t,batTotal,'--o');
title(strcat('ADMM Predication at Time',num2str(t0),' (Battery Capacity = ',num2str(mpc.B_upbound*Nbus),' p.u.)'));
xlabel('Time (Hour)');
ylabel('Power (p.u.)');
set(gca,'Xtick',t0:2:t0+T);
xlim([t0 t0+T]);
grid on;
legend('Sub Station Generation','Total Demand','Total Storage','Location','southeast');
figPredName = strcat('consensus_predication_time_',num2str(t0),'.png');
print(figPred,figPredName,'-dpng');

%% plot convergence 
figConv = figure;
xx = 1:length(objTotal);
subplot(2,1,1);
plot(xx,res.prim,'r-',xx,res.dual,'g-');
set(gca,'yscale','log');
legend('Primal Residual',...
    'Dual Residual','Location','southwest');
grid on;
xlim([1 length(objTotal)]);
ylim([1e-4 1]);
xlabel('Iteration');
ylabel('Residual');
title(['Consensus Predication at Time',num2str(t0),' Residual Convergence (\rho=' num2str(rho) ')']);

subplot(2,1,2);
yy = objTotal;
plot(xx,yy);
hold on;
plot([xx(1) xx(end)],[objTarget objTarget],':k','Linewidth',1);
set(gca,'yscale','log');
legend('ADMM Objective','Target');
grid on;
xlim([1 length(objTotal)]);
ylim([objTarget*(1-1e-2) objTarget*(1+1e-2)]);
xlabel('Iteration');
ylabel('Substation Power (p.u.)');
title(['Consensus Predication at Time',num2str(t0),' Objective Convergence (\rho=' num2str(rho) ')']);
set(figConv, 'Position', [100 100 1200 700]);
figstr = strcat('consensus_predication_convergence_time_',num2str(t0),'_rho_',num2str(rho));
%figname = strcat(figstr,'.fig');
%savefig(figname);
pngname = strcat(figstr,'.png');
print(figConv,pngname,'-dpng');
close all;

%% save calculation result
format shortG;
disp('--------------------------------');
disp(strcat('Prediction at T= ',num2str(t0),': Consensus Solved.'));
fprintf('Optimal value (Gurobi): %f \n', objTarget);
fprintf('Optimal value (Consensus ADMM): %f \n', objTotal(end));
fprintf('Exactness (Consensus ADMM): [%f %f]\n', x_check);


%% post infeasible conversion for early termination
%{

% closed-loop control input from distributed opf
s = complex(x_p,x_q);
s(1,:) = s(1,:) - complex(p_station,q_station);

% >>>> Version Adaption <<<< 
R_Line = mpc.branch(:,3)/mpc.base_Z;
X_Line = mpc.branch(:,4)/mpc.base_Z;

yalmip('clear');

clear S_sum;
[Nbus,T] = size(p_load_rated);
Nline = Nbus -1;

clear v S L s_station

% *********************** yalmip ************************************ 
% --- variables ---

s_station = sdpvar(1,T,'full','complex');

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
%objective = [];
objective = sum(sum(real(s_station).*real(s_station)));

% --- constraints ---
constraints = [];
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
        %constraints = constraints + [S(k,tt).^2 <= L(k,tt)*v(i,tt)];
    end
end
constraints = constraints + [ s_station + s(1,:) + S_sum(1,:) == 0];       % s_station + s(1) = net power at Bus 1

% voltage constraint    
constraints = constraints + [v(1,:) == 1];
constraints = constraints + [((1-mpc.dev_V))^2 <= v(2:end,:) <= ((1+mpc.dev_V))^2]; 


% --- options ---
ops = sdpsettings('verbose',0,'solver','gurobi');
%ops = sdpsettings('solver','gurobi','gurobi.FeasibilityTol',1e-9,'gurobi.OptimalityTol',1e-9,'verbose',1);
% --- solve --- 
%ctrloptmizer = optimizer(constraints, objective, ops, {b0,p_load_rated,q_load_rated,t0,bnow},{s_station,b,objective,L,v,s});

sol = optimize(constraints, objective, ops);
if sol.problem == 0
 % Extract and display value
 solution = value(s_station);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

p_station = real(double(s_station));
q_station = imag(double(s_station));

post_exact = abs(double(S)).^2 ./ (double(v(2:end,:)).*double(L));
%}

%% output 
sol.station0 = p_station(1,1);
%sol.b0 = batLocal(:,2);
sol.b = batLocal;
sol.x = x;
sol.z = z;
sol.lamda = lamda;

solfile = strcat('ADMM_sol_consensus_predication_time_',num2str(t0),'_rho_',num2str(rho),'.mat');

save(solfile,'sol','p_station','q_station','p_load_rated','q_load_rated','x_L','x_v','x_p','x_q','x_exact','batTotal','res','objTotal','objTarget','z');

end