function sol_opf = distr_opf_solver(net,sol_socp,p_pred,q_pred)

% :: closed-loop grid operation with distributed nonconvex algorithm 

for ii = 1:net.Ncalc

DoOnline = 1;    
    
t0 = ii-1;
sol_ref = sol_socp{ii};

net.p_demand = p_pred{ii};
net.q_demand = q_pred{ii};
if ii > 1
    net.e_bat_now = sol_opf{ii-1}.e_bat(:,1);
elseif ii == 1
    net.e_bat_now = net.e_bat_midnight;
end

% computation parameter
MIN_INN_EPS = 1e-4; 
MAX_RHO = 1e6;
RESCALE  = 1e-2;

if ii == 1
    MIN_EPS = 1e-4;
    if DoOnline == 0 
        RHO = 1;
        BETA = 1.1;
        MAX_ITER = 50;
        MAX_INN_ITER = 1000;
    elseif DoOnline == 1
        RHO = 1;
        BETA = 1.1;
        MAX_ITER = 100;
        MAX_INN_ITER = 1000;
    end
else
    MIN_EPS = 1e-4;
    RHO = 1;
    BETA = 1;
    MAX_ITER = 10;
    MAX_INN_ITER = 500;
end

INIT_TYPE = t0;

%% Initialization 
% parameter 
par.p_demand = net.p_demand;
par.q_demand = net.q_demand;

% variable 
switch INIT_TYPE
    case 0
        if DoOnline == 0
            [x,z,dual] = init_socp_warm(sol_ref,net);
        elseif DoOnline == 1
            [x,z,dual] = init_socp_warm(sol_ref,net);
            % [x,z,dual] = init_cold(sol_ref,net);            
        end
    otherwise
        [x,z,dual] = init_opf_warm(sol_opf{ii-1}.x,sol_opf{ii-1}.z,sol_opf{ii-1}.dual,net);
end

niter = 1;
clear resPrim resDual resPF
resPrim(niter,1) = 1;
resDual(niter,1) = 1;
resPF(niter,1) = 1;
obj = 0;

figure;
xx = 1:niter;
subplot(2,1,1);
hPrim = plot(xx,resPrim);
hold on;
hDual = plot(xx,resDual);
hold on;
hPF = plot(xx,resPF);
ylim([MIN_EPS 1]);
legend('primal','dual','PF');
ylabel('residual');
xlabel('niter');
set(gca,'yscale','log');
subplot(2,1,2);
hObj = plot(xx,obj);
hold on;
hRef = plot([1,xx(end)],[sol_ref.obj,sol_ref.obj],'k--');
legend('obj','socp ref');
% ylim([(sol_ref.obj*0.95) (sol_ref.obj*1.05)]);
ylabel('obj');
xlabel('niter');
drawnow;

%% parsing 

for idx = 1:net.Nbus
    funf{idx} = FunF(idx,net,RESCALE);
    gradf{idx} = GradF(idx,net,RESCALE);
    proxg{idx} = ProxG(idx,net);
end

%% Distributed non-convex opf solver 
while resPrim(end) >= MIN_EPS || resPF(end) >= MIN_EPS 
    
    % homotopy parameter 
    %{
    if ii > 1
        par.p_prev_pred = circshift(p_pred{ii-1},-1,2);
        par.p_demand = net.p_demand * niter / MAX_ITER + par.p_prev_pred * (MAX_ITER-niter) / MAX_ITER;
    end
    %}
    
    % primal-update 
    j = 0;
    x_res_Inn = 1;
    z_res_Inn = 1;
    z_old = z;
    while x_res_Inn > (MIN_INN_EPS/RHO) || z_res_Inn > (MIN_INN_EPS/RHO) || j < 3
        % x-update
        for idx = 1:net.Nbus 
            x_new{idx} = distr_x_opt(idx,x{idx},z,dual,net,par,RHO,t0,funf{idx},gradf{idx},proxg{idx});
        end
        % z-update
        for idx = 1:net.Nbus
            % z-optimization
            z_new{idx} = distr_z_opt(idx,x_new,z{idx},dual,net,RHO);
        end
        j = j + 1;
        [x_res_Inn,z_res_Inn] = resInnComp(x_new,z_new,x,z,net);
        x = x_new;
        z = z_new;
        if j >= MAX_INN_ITER
            fprintf('Inner loop truncated, niter = %d\n',niter);
            break;
        end
    end

    fprintf('j = %d\n',j);
    % dual-update
    for idx = 1:net.Nbus
        % dual-compute
        [dual_new{idx},res_sum_square(idx,1:2)] = distr_dual_update(idx,x_new{idx},z_new,dual{idx},net,par,RHO);
    end
    % residual
    resPF(niter,1) = sqrt(sum(res_sum_square(:,1)));
    resPrim(niter,1) = sqrt(sum(res_sum_square(:,2)));
    resDual(niter,1) = norm(cell2mat(z_new)-cell2mat(z_old),2);
    
    obj(niter,1) = 0;
    for idx = 1:net.Nbus
        if net.gen(idx,1) > 0 
            Nabus = net.abus(idx,1);
            if net.ncost == 3
                obj(niter,1) = obj(niter,1) + sum( net.c_gen(idx,1) * (x_new{idx}.bus(1+Nabus*3+1,:)).^2 );
            elseif net.ncost == 2
                obj(niter,1) = obj(niter,1) + sum( net.c_gen(idx,1) * (x_new{idx}.bus(1+Nabus*3+1,:) * net.baseMVA) );
            end
        end
    end
    
    % update
    x = x_new;
    z = z_new;
    dual = dual_new;
    
    if mod(niter,2) < 1 
        xx = 1:niter;
        set(hPrim,'XData',xx,'YData',resPrim);
        set(hDual,'XData',xx,'YData',resDual);
        set(hPF,'XData',xx,'YData',resPF);
        set(hObj,'XData',xx,'YData',obj);
        set(hRef,'XData',[1,niter]);
        refreshdata;
        drawnow;
    end
    
    if niter >= MAX_ITER 
        break;
    end 
    
    niter = niter + 1;
    
    if RHO <= MAX_RHO
        RHO = RHO * BETA;
    end
end

if niter >= MAX_ITER 
    disp('primal-dual not converge yet');
end    

xx = 1:length(resPrim);
set(hPrim,'XData',xx,'YData',resPrim);
set(hDual,'XData',xx,'YData',resDual);
set(hPF,'XData',xx,'YData',resPF);
set(hObj,'XData',xx,'YData',obj);
set(hRef,'XData',[1,niter]);
refreshdata;
drawnow;

% return;
%% Post Process  
for idx = 1:net.Nbus
    sol.v(idx,:) = double(x{idx}.bus(1,:));
    Nabus = net.abus(idx,1);
    Ngen = net.gen(idx,1);
    Nbat = net.bat(idx,1);
    if Nabus > 0 
        sol.L(idx,:) = double(x{idx}.bus(2,:));
        sol.P(idx,:) = double(x{idx}.bus(3,:));
        sol.Q(idx,:) = double(x{idx}.bus(4,:));
    end
    if Ngen > 0 
        sol.p_gen(idx,:) = double(x{idx}.bus(1+Nabus*3+1,:));
        sol.q_gen(idx,:) = double(x{idx}.bus(1+Nabus*3+2,:));
    end
    if Nbat > 0 
        sol.p_bat(idx,:) = double(x{idx}.bus(1+Nabus*3+Ngen*2+1,:));
        sol.e_bat(idx,:) = double(x{idx}.bus(1+Nabus*3+Ngen*2+2,:)); 
    end
end
sol.obj = obj(end);
sol.x = x;
sol.z = z;
sol.dual = dual;
sol.resPrim = resPrim;
sol.resDual = resDual;
sol.resPF = resPF;
sol.obj = obj;
sol.BETA = BETA;
format shortG


figure;
simplot(net,sol,t0);
hold on;
simplot(net,sol_ref,t0);

sol_opf{ii} = sol;

end

end