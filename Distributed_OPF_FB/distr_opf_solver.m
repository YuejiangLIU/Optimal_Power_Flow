%function [sol,x,z,dual] = distr_opf_solver(net)

% computation parameter
MAX_ITER = 600;
MIN_EPS = 1e-4;
PRIM_ITER = 20;

%% Model Parsing 
% :: min   fi(xi) + gi(xi) 
% ::  qp   fi(xi) = xi' * Fi * xi
% ::  lp   fi(xi) = Fi' * xi
% ::       g(xi) = Ic(xi), LBi <= xi <= UBi 
% :: s.t.  Ai * xi == bi 
% ::       xi' * Ei * xi == 0
% ::       xi == zi

Tpred = 1;
Cst.baseMVA = net.baseMVA;
for i = 1:net.Nbus
    % Constant
    Cst.cost = net.c_gen(i,:);
    Cst.ri = net.ri(i);
    Cst.xi = net.xi(i);
    Cst.bounds = [net.v_min(i),net.v_max(i);net.p_gen_min(i),net.p_gen_max(i);net.q_gen_min(i),net.q_gen_max(i)];
    Ncbus = net.cbus(i,1);
    if Ncbus > 0 
        for ic = 1:Ncbus
            Cst.rc(ic,1) = net.ri(net.cbus(i,ic+1));
            Cst.xc(ic,1) = net.xi(net.cbus(i,ic+1));
            if ic > 1
                error('please check ic ');
            end
        end
    end
    % matrix 
    [Ai,bi,Di,Ei,Fi,LBi,UBi] = pre_opf_matrix(net.gen(i),net.abus(i,1),net.cbus(i,1),Tpred,Cst);
    mat{i}.A = Ai;
    mat{i}.b = bi;
    mat{i}.D = Di;
    mat{i}.E = Ei;
    mat{i}.F = Fi;
    mat{i}.LB = LBi;
    mat{i}.UB = UBi;
end

%% Initialization 
rho = 1;

for i = 1:net.Nbus
    Nele = size(mat{i}.LB,1); 
    [x{i},z{i},dual{i}] = init_distri(Nele,net.gen(i),Tpred);    
end

for i = 1:net.Nbus
    % Z_hat 
    z_hat{i} = distr_commu(net.abus(i,:),net.cbus(i,:),z,i);
    % Parameter 
    Par.pdemand = net.p_demand(i);
    Par.qdemand = net.q_demand(i);
    % Matrix & Initialization Update
    mat{i}.b(2,1) = Par.pdemand;
    mat{i}.b(3,1) = Par.qdemand;
end

%% Distributed Computation
% format shortG
close all;
niter = 1;
clear resPrim resDual
resPrim(niter,1) = 1;
resDual(niter,1) = 1;
obj = 0;
% z_hat{i} = distr_commu(net.abus(i,:),net.cbus(i,:),z,i);

figure;
xx = 1:niter;
subplot(2,1,1);
hPrim = plot(xx,resPrim);
hold on;
hDual = plot(xx,resDual);
ylim([1e-4 1]);
ylabel('residual');
xlabel('niter');
set(gca,'yscale','log');
subplot(2,1,2);
hObj = plot(xx,obj);
ylabel('obj');
xlabel('niter');
drawnow

%% 
while resPrim(end) >= MIN_EPS || resDual(end) >= MIN_EPS
    % primal-update 
    j = 0;
    resInner = 1;
    z_old = z;
    while resInner > 1e-6 % for j = 1:PRIM_ITER
        % x-update
        for i = 1:net.Nbus %net.Nbus % net.Nbus
            % x-optimization
            x_new{i} = distr_x_opt(x{i},mat{i},z_hat{i},dual{i},rho);
        end
        % z-update
        for i = 1:net.Nbus
            % z-optimization
            z_new{i} = distr_z_opt(x_new,net.abus,net.cbus,dual,i,net.gen,rho);
        end
        for i = 1:net.Nbus
            % z-communication
            z_hat{i} = distr_commu(net.abus(i,:),net.cbus(i,:),z_new,i);
        end
        resInner = max(vec(cell2mat(z_new))./vec(cell2mat(z_old)));
        z_old = z_new;
        j = j + 1;
        if j >= PRIM_ITER
            break;
        end
    end
    fprintf('j = %d\n',j);
    % dual-update
    for i = 1:net.Nbus
        % dual-compute
        dual_new{i} = distr_dual_update(x_new{i},mat{i},z_hat{i},dual{i},rho);
    end
    % residual
    resPrimSumSqu = 0;
    for i = 1:net.Nbus
        resPrimSumSqu = resPrimSumSqu + (mat{i}.D * x_new{i} - z_hat{i})'*(mat{i}.D * x_new{i} - z_hat{i});
    end
    resPrim(niter,1) = sqrt(resPrimSumSqu);
    resDual(niter,1) = rho * norm(cell2mat(z_new)-cell2mat(z),2);
    obj(niter,1) = (mat{1}.F' * x{1} + mat{2}.F' * x{2}) * net.baseMVA;
    % update
    x = x_new;
    z = z_new;
    dual = dual_new;

    if mod(niter,10) < 1
        xx = 1:niter;
        set(hPrim,'XData',xx,'YData',resPrim);
        set(hDual,'XData',xx,'YData',resDual);
        set(hObj,'XData',xx,'YData',obj);
        refreshdata
        drawnow
    end
    
    if niter >= MAX_ITER
        break;
    end    
    niter = niter + 1;
    
end

if niter >= MAX_ITER 
    disp('primal-dual not converge yet');
end    

%% Post Process  
for i = 1:net.Nbus
    sol.v(i,1) = double(x{i}(1));
    sol.L(i,1) = double(x{i}(2));
    sol.P(i,1) = double(x{i}(3));
    sol.Q(i,1) = double(x{i}(4));
    if net.gen(i)
        sol.p(i,1) = double(x{i}(5)) * net.gen(i);
        sol.q(i,1) = double(x{i}(6)) * net.gen(i);
    end
end
sol.obj = obj(end);
    
