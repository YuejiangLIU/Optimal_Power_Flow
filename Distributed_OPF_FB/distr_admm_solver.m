%function distr_admm_solver(net)

MAX_ITER = 200;
RHO = 10;

yalmip('clear');

% Parameter
p_load_rated = - net.p_demand;
q_load_rated = - net.q_demand;

T = 1;

clear x z lamda

%% Initialization 
for ii = 1:net.Nbus
    % --- local itself ---
    x(ii).bus = ones(6,1);
    lamda(ii).bus = ones(4,1);
    % --- local copy of ancestor --
    Nabus = net.abus(ii,1);
    if Nabus>0
        x(ii).abus = ones(1,1);              % [ ind of A; v_A];
        lamda(ii).abus = ones(1,1);
    end
    % --- local copy of children ---
    Ncbus = net.cbus(ii,1);
    if Ncbus>0 
        x(ii).cbus = ones(3,Ncbus);             % [ ind of C; L_C, P_C, Q_C];
        lamda(ii).cbus = ones(3,Ncbus);
    end     
end
% z initialization
z = ones(4,net.Nbus);

%%
%figure;
Ncurrentiter = 0;
for niter = Ncurrentiter + 1:Ncurrentiter + MAX_ITER
    %% x-update
    for ii = 1:net.Nbus
        Nabus = net.abus(ii,1);
        Ncbus = net.cbus(ii,1);
        if niter == 1
            yalmip('clear');
            clear lamdaAA lamdaBB lamdaCC ZZ
            lamdaBB = sdpvar(4,1,'full');
            if Nabus > 0
                lamdaAA = sdpvar(1,Nabus,'full');
            else
                clear lamdaAA;
            end
            if Ncbus > 0
                lamdaCC = sdpvar(3,Ncbus,'full');
            else
                clear lamdaCC;
            end
            ZZ = sdpvar(4,net.Nbus,'full');
            
            % ---- variable ----
            clear xB;
            xB = sdpvar(6,1,'full');    % [v; L; P; Q; p; q]
            
            clear xA;
            if Nabus>0
                xA = sdpvar(1,1,'full');             % [va]
            else
                xA = [];
            end
            
            clear xC;
            if Ncbus>0
                xC = sdpvar(3,Ncbus,'full');       % [Lc; Pc; Qc]
            else
                xC = [];
            end
            
            %BAT = sdpvar(1,T+1,'full');
            %p_battery = sdpvar(1,T,'full');
            
            % ---- auxiliary variable ----
            p_generator = sdpvar(1,1,'full');
            q_generator = sdpvar(1,1,'full');
            % q_capacitor = sdpvar(1,T,'full');
            
            
            % ---- expression in function and constraint ----
            % >>> sub objective <<<
            clear fi;
            % fi = norm(p_generator,2);               % flatten generation variance
            fi = net.c_gen(ii,1) * p_generator;
            
            clear productB;
            productB = sum(sum(lamdaBB.*xB(1:4,1)));
            
            clear productA;
            if Nabus > 0
                productA = sum(lamdaAA.*xA);
            else
                productA = 0;
            end
            
            clear productC;
            if Ncbus > 0
                if Ncbus == 1
                    productC = sum(sum(lamdaCC.*xC));
                else        % Nc > 1 => xC is a matrix
                    productC = sum(sum(sum(lamdaCC.*xC)));
                end
            else
                productC = 0;
            end
            
            clear sumsqB;
            sumsqB = sum(sum_square(xB(1:4,1) - ZZ(:,ii)));
            
            clear sumsqA;
            if Nabus > 0
                indA = net.abus(ii,2);
                sumsqA = sum_square( xA - ZZ(1,indA) );
            else
                sumsqA = 0;
            end
            
            clear sumsqC;
            sumsqC = 0;
            if Ncbus > 0
                for cc = 1:Ncbus               % enumerate chilren
                    indC = net.cbus(ii,cc+1);      % who is the child
                    sumsqC = sumsqC + sum(sum_square( xC(:,cc) - ZZ(2:4,indC)));
                end
            else
                sumsqC = 0;
            end
            
            clear sumPC sumQC;
            sumPC = 0;
            sumQC = 0;
            if Ncbus > 0
                for cc = 1:Ncbus
                    indC = net.cbus(ii,cc+1);
                    sumPC = sumPC + ( xC(2,cc) - xC(1,cc) * net.ri(indC) );
                    sumQC = sumQC + ( xC(3,cc) - xC(1,cc) * net.xi(indC) );
                end
            else
                sumPC = 0;
                sumQC = 0;
            end
            
            % ---- objective ----
            objective = fi + productB + productA + productC ...
                + RHO/2 * ( sumsqB + sumsqA + sumsqC );
            
            % ---- constraint ----
            constraints = [];
            % power flow constraints
            if Nabus > 0 
                constraints = constraints + [xA == xB(1,1) - 2*(net.ri(ii)*xB(3,1)+net.xi(ii)*xB(4,1)) + xB(2,1)*(net.ri(ii)^2+net.xi(ii)^2) ];
            else
                % constraints = constraints + [xB(1,1) == 1];
                constraints = constraints + [xB(2,1) == 0];
                constraints = constraints + [xB(3,1) == 0];
                constraints = constraints + [xB(4,1) == 0];
            end
            constraints = constraints + [sumPC + xB(5,1) == xB(3,1)];
            constraints = constraints + [sumQC + xB(6,1) == xB(4,1)];
            % voltage constraints
            constraints = constraints + [net.v_min(ii)^2 <= xB(1,1) <= net.v_max(ii)^2];
            % power injection constraints
            constraints = constraints + [xB(5,1) == p_load_rated(ii) + p_generator ];
            constraints = constraints + [xB(6,1) == q_load_rated(ii) + q_generator ];
            if net.gen(ii) > 0
                constraints = constraints + [net.p_gen_min(ii) <= p_generator <= net.p_gen_max(ii)];
                constraints = constraints + [net.q_gen_min(ii) <= q_generator <= net.q_gen_max(ii)];
            else
                constraints = constraints + [p_generator == 0];
                constraints = constraints + [q_generator == 0];
            end
            
            % ---- settings ----
            ops = sdpsettings('solver','gurobi');
            
            % ---- modeling ----
            if Nabus < 1
                xoptimizer{ii} = optimizer(constraints, objective, ops,...
                    {lamdaBB,lamdaCC,ZZ},{xB,xC,p_generator,q_generator});
            elseif Ncbus < 1
                xoptimizer{ii} = optimizer(constraints, objective, ops,...
                    {lamdaAA,lamdaBB,ZZ},{xA,xB,p_generator,q_generator});
            else
                xoptimizer{ii} = optimizer(constraints, objective, ops,...
                    {lamdaAA,lamdaBB,lamdaCC,ZZ},{xA,xB,xC,p_generator,q_generator});
            end
        end
        
        % --- YALMIP Controller ---
        if Nabus < 1
            xctrl = xoptimizer{ii};
            [xopt, errcode] = xctrl{{lamda(ii).bus,lamda(ii).cbus,z}};
            xB = xopt{1};
            xC = xopt{2};
            p_generator = xopt{3};
            q_generator = xopt{4};
        elseif Ncbus <1
            xctrl = xoptimizer{ii};
            [xopt, errcode] = xctrl{{lamda(ii).abus,lamda(ii).bus,z}};
            xA = xopt{1};
            xB = xopt{2};
            p_generator = xopt{3};
            q_generator = xopt{4};
        else
            xctrl = xoptimizer{ii};
            [xopt, errcode] = xctrl{{lamda(ii).abus,lamda(ii).bus,lamda(ii).cbus,z}};
            xA = xopt{1};
            xB = xopt{2};
            xC = xopt{3};
            p_generator = xopt{4};
            q_generator = xopt{5};
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
        
        x(ii).bus = double(xB);
        if Nabus>0
            x(ii).abus = double(xA);
        end
        if Ncbus>0
            x(ii).cbus = double(xC);
        end
        
        p_gen(ii,:) = double(p_generator);
        q_gen(ii,:) = double(q_generator);
        objLocal(ii,1) = net.c_gen(ii,1) * double(p_generator);
        
    end
    
    objTotal(niter,1) = sum(objLocal) * net.baseMVA;

    %% ADMM z-update
    zold = z;
    for ii = 1:net.Nbus
        
        Nabus = net.abus(ii,1);
        Ncbus = net.cbus(ii,1);
        
        % --- update z_v ---
        if Ncbus > 0
            facA = RHO/2 * (1 + Ncbus);
            facB = - RHO*x(ii).bus(1) - lamda(ii).bus(1);
            for cc = 1:Ncbus
                indC = net.cbus(ii,cc+1);
                facB = facB - RHO*x(indC).abus(1) - lamda(indC).abus(1);
            end
            z(1,ii) = -facB/facA/2;
        else
            z(1,ii) = x(ii).bus(1) + lamda(ii).bus(1)/RHO;
        end
        clear facA;
        clear facB;
        
        % --- update z_L, z_P, z_Q ---
        if Nabus > 0
            indA = net.abus(ii,2);         % who is ii's ancestor
            nbC = find(net.cbus(indA,:) == ii,1,'last') - 1; % ii is his ancestor's No. # child
            facA = RHO;
            facB(1,1) = RHO*(-x(ii).bus(2,1) - x(indA).cbus(1,nbC))-lamda(ii).bus(2,1) - lamda(indA).cbus(1,nbC);
            facB(2,1) = RHO*(-x(ii).bus(3,1) - x(indA).cbus(2,nbC))-lamda(ii).bus(3,1) - lamda(indA).cbus(2,nbC);
            facB(3,1) = RHO*(-x(ii).bus(4,1) - x(indA).cbus(3,nbC))-lamda(ii).bus(4,1) - lamda(indA).cbus(3,nbC);
            z(2,ii) = -facB(1,1)/facA/2;
            z(3,ii) = -facB(2,1)/facA/2;
            z(4,ii) = -facB(3,1)/facA/2;
        else
            z(2:4,ii) = x(ii).bus(2:4,1) + lamda(ii).bus(2:4,1)/RHO;
        end
        clear facA;
        clear facB;
    end
    
    %% Multiplier update  
    for ii = 1:net.Nbus  
        Nabus = net.abus(ii,1);
        Ncbus = net.cbus(ii,1);
        
        lamda(ii).bus = lamda(ii).bus + RHO * ( x(ii).bus(1:4,1) - z(:,ii) );
        
        if Nabus > 0
            indA = net.abus(ii,2);
            lamda(ii).abus = lamda(ii).abus + RHO * ( x(ii).abus - z(1,indA) );
        end
        
        if Ncbus > 0
            for cc = 1:Ncbus
                indC = net.cbus(ii,cc+1);
                lamda(ii).cbus(:,cc) = lamda(ii).cbus(:,cc) + RHO * ( x(ii).cbus(:,cc) - z(2:4,indC) );
            end
        end       
    end
    
    %% residual
    res.bus(niter,1) = 0;
    for ii = 1:net.Nbus
        %res.B(niter,1) = res.B(niter,1) + norm( x(ii).B-z(:,ii,:),2 );
        res.bus(niter,1) = res.bus(niter,1) + sum(sum_square(x(ii).bus(1:4,1)-z(:,ii)));
    end
    
    res.abus(niter,1) = 0;
    res.cbus(niter,1) = 0;
    for ii = 1:net.Nbus
        Nabus = net.abus(ii,1);
        if Nabus > 0
            indA = net.abus(ii,2);
            res.abus(niter,1) = res.abus(niter,1) + sum_square( x(ii).abus - z(1,indA));
        end
        
        Ncbus = net.cbus(ii,1);
        if Ncbus > 0
            for cc = 1:Ncbus
                indC = net.cbus(ii,cc+1);
                res.cbus(niter,1) = res.cbus(niter,1) + sum(sum_square( x(ii).cbus(:,cc) - z(2:4,indC)));
            end
        end
    end
    
    resPrimal = sqrt(res.abus(niter) + res.bus(niter) + res.cbus(niter));
    res.prim(niter,1) = resPrimal;
    resDual = RHO*sqrt(sum(sum(sum_square(z-zold))));
    res.dual(niter,1) = resDual;
    
end

xx = 1:niter;
figure;
subplot(2,1,1);
plot(res.prim);
hold on;
plot(res.dual);
legend('primal','dual');
ylabel('residual');
xlabel('niter');
set(gca,'yscale','log');
subplot(2,1,2);
plot(objTotal);
ylabel('obj');
xlabel('niter');
    
Ncurrentiter = niter;
for ii = 1:net.Nbus    
x_v(ii,1) = x(ii).bus(1,1);
x_L(ii,1) = x(ii).bus(2,1);
x_P(ii,1) = x(ii).bus(3,1);
x_Q(ii,1) = x(ii).bus(4,1);
x_p(ii,1) = x(ii).bus(5,1);
x_q(ii,1) = x(ii).bus(6,1);
end

z_v = z(1,:);
z_L = z(2,:);
z_P = z(3,:);
z_Q = z(4,:);