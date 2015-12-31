%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function distopf(PERTUNIT,Niter,RHO) 

clear all; close all; clc; 

gadd = 0;
dev = [0.1 0.1];
mode = 1;

%% parameter setup

PERTUNIT = 4;           % perturbation percentage
Niter = 100;            % max iteration of distributed optimization 
RHO = 1;              % penalty parameter of distributed optimization

Npred = 1;             % number of sampling times in closed-loop simulation
plotind = [];           % list of plots for demand prediction
scalefactor = 1;        % scale factor of generator & capacitor (ref: Gan, L., etc (2015). Exact Convex Relaxation of Optimal Power Flow in Radial Networks. IEEE Transactions on Automatic Control 60, 72–87.)

battery_capacity = 0;        % battery capacity assumed at each bus 
battery_initial = 1/5;          % battery charged percentage at each bus at time 0  

% flag_result = 1;                                

DoCentComp = 1;     % Do centralized optimization or not? (1 = Yes, 0 = No)

HIDEFIG = 0;        % Hide figures or not during program execution? (1 = Yes, 0 = No) 
if HIDEFIG
    set(0,'DefaultFigureVisible','off');        
end

% range of perturbation
PERTMIN = PERTUNIT;         
PERTMAX = PERTUNIT;        
pertmin = PERTMIN/100;            
pertmax = PERTMAX/100;

% type of distributed algorithm employed  
TypeAlg = 'BCD';        
%TypeAlg = 'ADMM'; 
%TypeAlg = 'GS'; 
%TypeAlg = 'Homotopy'; 

%% file and folder
%{
clear foldername;
rootDir = pwd;                  % directory of root folder 
foldername = strcat(TypeAlg,'_Pert',num2str(PERTUNIT),'_Iter',num2str(Niter),'_Rho',num2str(RHO));      % name of new folder

if exist(foldername,'dir') > 0
    rmdir(foldername,'s');      % delete folder >> Attention! <<
end 

mkdir(foldername);              % create new folder for new simulation 

folderDir = strcat(rootDir,'\',foldername); % directory of sub folder ( '\' might be different on different OS)

% copy scripts to sub folder 
copyfile('distopf.m',folderDir);
copyfile('loadpred.m',folderDir);
copyfile('network56.m',folderDir);
copyfile('network7.m',folderDir);
copyfile('centralized_controller.m',folderDir);
%copyfile('sol_dayahead_7bus_rho_0.1.mat',folderDir);
filecontroller = strcat('distributed_',TypeAlg,'_controller.m'); 
copyfile(filecontroller,folderDir);

centfile = strcat('centralized_mpc_sol_pert_',num2str(PERTUNIT),'.mat');
if exist(centfile)
    copyfile(centfile,folderDir);
end

if exist('sol_consensus_predication_time_offline.mat')
    copyfile('sol_consensus_predication_time_offline.mat',folderDir);
end

cd(folderDir);

disp('==============================');
disp(foldername);           % display folder / simulation setup
%}

%% raw data 

% read newtork data 
% mpc = network7;             
%mpc = network56(gadd,dev);    

mpc = case14Tree;



return;

% data pre-processing
T = 1;
if mode > 1
demand = repmat(mpc.injection(:,2),1,T).*mpc.demand(:,2:T+1);
elseif mode == 1
demand = mpc.injection(:,2);
end
generator = mpc.injection(:,3);
capacitor = mpc.injection(:,4);
battery = mpc.injection(:,5);

%% demand prediction and perturbation 

disp('---------------------------------------');
disp('Status: Online load is to predict ... ');

demandpred = loadpred(demand,Npred,plotind,pertmin,pertmax);


%% network info processing  

% topology 
From_Bus = mpc.branch(:,1);
To_Bus = mpc.branch(:,2);
R_Line = mpc.branch(:,3)/mpc.base_Z;

% total number of buses/lines
mpc.Nbus = length(mpc.injection);       
mpc.Nline = length(mpc.branch);         

% impedance
R = zeros(mpc.Nbus,1);
X = zeros(mpc.Nbus,1);
Z = complex(R,X);
for kk = 1:mpc.Nline
    i = mpc.branch(kk,1);
    h = mpc.branch(kk,2);
    R(h) = mpc.branch(kk,3);
    X(h) = mpc.branch(kk,4);
    Z(h) = complex(R(h),X(h));
end
mpc.R = R./mpc.base_Z;
mpc.X = X./mpc.base_Z;
mpc.Z = Z./mpc.base_Z;

% constraints on device power injection 
mpc.p_generator_lowbound = repmat(zeros(mpc.Nbus,1),1,T);
mpc.s_generator_upbound = repmat(generator*scalefactor,1,T);

mpc.p_capacitor_rated = repmat(zeros(mpc.Nbus,1),1,T);
mpc.q_capacitor_upbound = repmat(capacitor*scalefactor,1,T);

mpc.B_upbound = battery_capacity * repmat(battery,1,T);

mpc.A = zeros(mpc.Nbus,1);
mpc.C = zeros(mpc.Nbus,1);

for k = 1:mpc.Nline
    j = mpc.branch(k,1);
    i = mpc.branch(k,2);
    mpc.C(j,1) = mpc.C(j,1) + 1;        % number of ancestors 
    mpc.C(j,mpc.C(j,1)+1) = i;          % index of ancestors 
    mpc.A(i,1) = mpc.A(i,1) + 1;
    mpc.A(i,mpc.A(i,1)+1) = j;
end


%% Centralized Optimization
if DoCentComp 
    disp('-------------------------------------------');
    disp('Status: Centralized optimization starts ... ');

    tStart=tic;

	% --- Centralized Offline --- 
	% offline setup 
	b0 = battery_initial*battery_capacity*battery;         % battery state at each bus at time 0 
	t0 = 0;
	bnow = b0;                                                      % battery state for the current scheduling prob 
	p_load_rated = - demand*mpc.pf;
	q_load_rated = - demand*sqrt(1-mpc.pf^2);
	% offline computation 
	ctrl = centralized_controller(mpc,b0,p_load_rated,q_load_rated,t0,bnow,mode);        % centralized offline 
	if mode > 1
        bOffline = ctrl.b;
    end
	objOffline = ctrl.obj;
	LOffline = ctrl.L;
	genOffline = ctrl.gen;
	lossOffline = sum(LOffline.*repmat(R_Line,1,T));
    exactOffline = ctrl.exact;
    if mode > 1
	% offline plot
	figOffline = figure;
	t = 0:1:T;
	plot(t,double([real(genOffline) real(genOffline(1))]),'--d',...
		t,-[sum(-demand*mpc.pf) sum(-demand(:,1)*mpc.pf)],'--*',...
		t,sum(bOffline),'--o');
	title(strcat('YALMIP Offline Prediction (Battery Capacity = ',num2str(sum(sum(mpc.B_upbound))),' p.u.)'));
	xlabel('Time (Hour)');
	ylabel('Power (p.u.)');
	grid on;
	set(gca,'Xtick',0:2:T);
	xlim([0 T]);
	%ylim([0 4]);
	legend('Generation','Total Demand','Total Storage','Location','southeast');
	figstr = strcat('sim_offline');
	pngname = strcat(figstr,'.png');
	print(figOffline,pngname,'-dpng');
	% close all;
    end
    exactMax = max(max((exactOffline)))
    
    max(max(ctrl.v./(1+dev(2)).^2))
    min(min(ctrl.v./(1-dev(1)).^2))
    return;    
    
    %%
	% --- Centralized Online ---  
	%plotind = 1:Npred;             % online plot index 
	for ii = 1:Npred  
		% online setup
		t0 = ii-1;
		p_load_rated = - demandpred{ii}*mpc.pf;
		q_load_rated = - demandpred{ii}*sqrt(1-mpc.pf^2);
		% online computation    
		ctrl = centralized_controller(mpc,b0,p_load_rated,q_load_rated,t0,bnow);
		b = ctrl.b;
		objective = ctrl.obj;
		L = ctrl.L;
		s_gen = ctrl.gen;
		v = ctrl.v;
		s = ctrl.s;
		% online plot
		if (any(ii == plotind))
			t = ii-1:1:T+ii-1;
			figSim = figure;
			plot(t,[real(s_gen) real(s_gen(1))],'--d',...
			t,-[sum(p_load_rated) sum(p_load_rated(:,1))],'--*',...
			t,sum(b),'--o');
			title(strcat('YALMIP Prediction at Time ',num2str(ii-1),' (Battery Capacity = ',num2str(battery_capacity*mpc.Nbus),' p.u.)'));
			xlabel('Time (Hour)');
			ylabel('Power (p.u.)');
			grid on;
			set(gca,'Xtick',ii-1:2:T+ii-1);
			xlim([ii-1 T+ii-1]);
			ylim([0 0.8]);
			legend('Sub Station Generation','Total Demand','Total Storage','Location','southeast');
			figstr = strcat('sim_time_',num2str(ii));
			%figname = strcat(figstr,'.fig');
			%savefig(figname);
			pngname = strcat(figstr,'.png');
			saveas(figSim,pngname,'png');
			close all;
		end
		% online results store
		b0Online(:,ii) = bnow;
		genOnline(1,ii) = s_gen(1);
		objOnline(ii,1) = objective;
		demandOnline(:,ii) = demandpred{ii}(:,1);
		lossOnline(1,ii) = sum(L(:,1).*R_Line);
		LOnline(:,ii) = L(:,1);
		vOnline(:,ii) = v(:,1);
		sOnline(:,ii) = s(:,1);
		% battery state update
		bnow = double(b(:,2));
	end
	
	% mpc siumulation 
	%{
	figMPC = figure;
	t = 0:1:T;
	plot(t,real(genOnline),'--d',...
		t,sum(demandOnline*mpc.pf),'--*',...
		t,sum(b0Online),'--o');
	hold on;
	plot([0 T],[battery_capacity*mpc.Nbus battery_capacity*mpc.Nbus],':k','Linewidth',1);
	hold off;
	title(strcat('YALMIP MPC Simulation (Battery Capacity = ',num2str(battery_capacity*mpc.Nbus),' p.u.)'));
	xlabel('Time (Hour)');
	ylabel('Power (p.u.)');
	grid on;
	set(gca,'Xtick',0:2:T);
	xlim([0 T]);
	%ylim([0 0.8]);
	legend('Sub Station Generation','Total Demand','Total Storage','Location','southeast');
	figstr = 'sim_MPC';
	pngname = strcat(figstr,'.png');
	saveas(figMPC,pngname,'png');
	close all;
	%}
	
	disp('-----------------------------------------------');
	disp('Status: Centralized optimization has been done. ');
	disp('------------------------------------------------');
	disp('Ready for the distributed optimizaiton now. ');
	% disp('Please change the DoCentComp to 0, and run the script agian.');
	
	ComputeTimeOptimizer=toc(tStart)
	
	save(centfile);
	
	cd(rootDir);
	save(centfile);    
% return;
    cd(folderDir);
else 
    disp('-----------------------------------------------------');    
    disp('Status: Centralized optimization has been loaded ... ');
    load(centfile,'-regexp', '^(?![RHO,TypeAlg])\w');
    b0 = battery_capacity/3*ones(mpc.Nbus,1);
end


%% Centralized Offline Simulation 

%{
genOpenLoop = [genOffline genOffline(1)];           
demandCloseLoop = demandOnline;                      
lossOpenLoop = [lossOffline lossOffline(1)];
powergapOpenLoop = real(genOpenLoop) - sum(demandCloseLoop*mpc.pf) - lossOpenLoop; 
bOpenLoop(1,1) = sum(b0);
for t=1:T
    sumgapOpenLoop(1,t) = sum(powergapOpenLoop(1:t));
end
bOpenLoop(1,2:T+1) = sum(b0) + sumgapOpenLoop(1:T);
% open loop simulation  
figOpen = figure;
t = 0:1:T;
plot(t,real(genOpenLoop),'--d',...
    t,sum(demandCloseLoop*mpc.pf),'--*',...
    t,bOpenLoop,'--o');
hold on;
plot([0 T],[battery_capacity*mpc.Nbus battery_capacity*mpc.Nbus],':k','Linewidth',1);
hold off;
title(strcat('Centralized Offline Simulation (Battery Capacity = ',num2str(battery_capacity*mpc.Nbus),' p.u.)'));
xlabel('Time (Hour)');
ylabel('Power (p.u.)');
grid on;
set(gca,'Xtick',0:2:T);
xlim([0 T]);
%ylim([0 0.8]);
legend('Sub Station Generation','Total Demand','Total Storage','Location','southeast');
figstr = 'sim_open_loop';
pngname = strcat(figstr,'.png');
print(figOpen,pngname,'-dpng');
close all;
%}

%% distributed optmization parameter 

DoOfflineDist = 0;                               % Do offline distributed optimization or not? (1 = Yes, 0 = No)

%% Distributed Algorithm Initialization (notation refers to my master thesis)
clear mpc.x mpc.lamda mpc.z
% x & lamda initialization
for ii = 1:mpc.Nbus
    % --- local itself ---
    mpc.x(ii).B = ones(6,1,T);
    mpc.lamda(ii).B = ones(4,1,T);
    % --- local copy of ancestor --
    Na = mpc.A(ii,1);
    if Na>0
        mpc.x(ii).A = ones(1,1,T);              % [ ind of A; v_A];
        mpc.lamda(ii).A = ones(1,1,T);
    end
    % --- local copy of children ---
    Nc = mpc.C(ii,1);
    if Nc>0 
        mpc.x(ii).C = ones(3,Nc,T);             % [ ind of C; L_C, P_C, Q_C];
        mpc.lamda(ii).C = ones(3,Nc,T);
    end     
end
% z initialization
mpc.z = ones(4,mpc.Nbus,T);

%% --- Offline Distributed Optimization ---
clear sol;
if ~DoOfflineDist 
    % Load offline data, if it has been calculated
    %solfile = strcat('sol_dayahead_7bus_rho_0.1.mat');
    solfile = strcat('sol_consensus_predication_time_offline.mat');
    load(solfile,'-regexp', '^(?![Do,RHO,TypeAlg])\w');
    disp('------------------------------------------');    
    disp('Status: Offline distributed opf data has been loaded... ');
else 
    solfile = strcat('sol_dayahead_rho_',num2str(RHO),'.mat');
    disp('----------------------------------');
    disp('Status: Offline distributed opf is running...');
    
    % offline computation
    t0 = 0;                
    bnow = b0;
    NiterOffline = 1000;
    sol = distributed_ADMM_controller(demand,b0,t0,bnow,mpc,RHO,NiterOffline,objOffline);
    save(solfile,'sol');

    %rename offline result
    loadfile = strcat('ADMM_sol_consensus_predication_time_',num2str(t0),'_rho_',num2str(RHO),'.mat');
    load(loadfile);

    offlinefile = strcat('sol_consensus_predication_time_offline.mat');
    save(offlinefile,'sol','res','objTotal','objTarget','p_station','p_load_rated','batTotal');

    disp('-----------------------------------');
    disp('Status: Offline distributed opf has been done.');
    disp('-----------------------------------');
	disp('Ready for the distributed optimizaiton now. ');
	% disp('Please change the DoCentComp to 0, and run the script agian.');
    % return; 
end

%% --- Online Distributed MPC ---   
disp('Status: Online Distributed OPF is to start.');
tStart=tic;
%Warm Starting
bnow = b0;
mpc.x = sol.x;
mpc.z = sol.z;
mpc.lamda = sol.lamda;
clear sol;
for t0 = 0:Npred-1
    % Online Problem
    demandnow = demandpred{t0+1};          % perturbance
    % Online Optimization 
    if TypeAlg == 'ADMM'
        sol = distributed_ADMM_controller(demandnow,b0,t0,bnow,mpc,RHO,Niter,objOnline(t0+1));  %return;
    elseif TypeAlg == 'BCD'
        sol = distributed_BCD_controller(demandnow,b0,t0,bnow,mpc,RHO,Niter,objOnline(t0+1));
    end
    return;
    
    % Battery update for next sampling time
    bnow = sol.b(:,2);
    % Warm Starting for next sampling time
    % --- warm x & lamda ---
    for ii = 1:mpc.Nbus
        mpc.x(ii).B(:,:,1:T-1) = sol.x(ii).B(:,:,2:T);
        mpc.lamda(ii).B(:,:,1:T-1) = sol.lamda(ii).B(:,:,2:T);
        mpc.x(ii).B(:,:,T) = sol.x(ii).B(:,:,1);
        mpc.lamda(ii).B(:,:,T) = sol.lamda(ii).B(:,:,1);
        
        Na = mpc.A(ii,1);
        if Na>0        
            mpc.x(ii).A(:,:,1:T-1) = sol.x(ii).A(:,:,2:T);
            mpc.lamda(ii).A(:,:,1:T-1) = sol.lamda(ii).A(:,:,2:T);
            mpc.x(ii).A(:,:,T) = sol.x(ii).A(:,:,1);
            mpc.lamda(ii).A(:,:,T) = sol.lamda(ii).A(:,:,1);
        end

        Nc = mpc.C(ii,1);
        if Nc>0 
            mpc.x(ii).C(:,:,1:T-1) = sol.x(ii).C(:,:,2:T);
            mpc.lamda(ii).C(:,:,1:T-1) = sol.lamda(ii).C(:,:,2:T);            
            mpc.x(ii).C(:,:,T) = sol.x(ii).C(:,:,1); 
            mpc.lamda(ii).C(:,:,T) = sol.lamda(ii).C(:,:,1);             
        end
    end
    % --- warm z ---
    mpc.z(:,:,1:T-1) = sol.z(:,:,2:T);
    mpc.z(:,:,T) = sol.z(:,:,1);  
    
end

ComputeTimeDist=toc(tStart)

%% Online Distributed OPF Analysis 
disp('--------------------------------');
disp('Status: Online Distributed OPF is to start.');

% read distributed opf result 
Npred = 25;
for ii = 1:Npred
    tt = ii - 1;
    filename = strcat(TypeAlg,'_sol_consensus_predication_time_',num2str(tt),'_rho_',num2str(RHO),'.mat');
    load(filename,'sol','p_station','p_load_rated','x_L','x_v','x_p','x_q','x_exact','batTotal','res','objTotal','objTarget');
    exact_dist(:,ii) = x_exact(:,1);
    objRatio(ii,1) = objTotal(end)/objTarget;
    primalStart(ii,1) = res.prim(1);
    primalEnd(ii,1) =  res.prim(end);
    dualStart(ii,1) = res.dual(1);
    dualEnd(ii,1) = res.dual(end);
    %p_station_dist(ii,1) = (x_p(1,1) - p_load_rated(1,1))+(x_p(45,1) - p_load_rated(45,1));
    %p_station_dist(ii,1) = sum(x_p(:,1) - p_load_rated(:,1));
    %p_gen_dist(ii,1) = p_gen(1);
    p_gen_dist(ii,1) = p_station(1);
    p_demand_dist(ii,1) = -sum(p_load_rated(:,1));
    bnow_dist(ii,1) = batTotal(1);
    L_dist(:,ii) = x_L(:,1);
    v_dist(:,ii) = x_v(:,1);
    p_dist(:,ii) = x_p(:,1);
    q_dist(:,ii) = x_q(:,1); 
end

bnowRatio = bnow_dist./sum(b0Online)';          % ratio of battery scheduling result from dist sol and optimal sol 
genRatio = p_gen_dist./real(genOnline)';        % ratio of geneation scheduling result from dist sol and optimal sol 

%{
% plot validatoin
T = Npred-1;
t = 0:T;
figValid = figure; 
% exactness 
subplot(4,1,1);
plot(t,max(exact_dist),'-o');
ylabel('Exact');
grid on;
xlim([0 Npred-1]);
ylim([1 1.1]);
% objective
subplot(4,1,2);
plot(t,objRatio,'-p',t,genRatio,'-*',t,bnowRatio,'-s');
legend('Objective Comparison','Generation Comparison','Storage Comparison','Location','west');
ylabel('Relative Ratio');
grid on;
xlim([0 Npred-1]);
ylim([1-1e-1 1+1e-1]);
% dual residual
subplot(4,1,3);
plot(t,primalStart,'-^',t,primalEnd,'-v');
legend('Primal Redisual Start','Primal Residual End','Location','west');
ylabel('Primal Residual');
grid on;
xlim([0 Npred-1]);
ylim([1e-4 1e-1]);
set(gca,'yscale','log');
% dual residual
subplot(4,1,4);
plot(t,dualStart,'-^',t,dualEnd,'-v');
legend('Dual Residual Start','Dual Residual End','Location','west');
ylabel('Dual Residual');
grid on;
xlim([0 Npred-1]);
ylim([1e-5 1e-2]);
set(gca,'yscale','log');
set(figValid, 'Position', [0 0 1200 700]);
print(figValid,'validation_consensus.png','-dpng');
%close(figExact);

% plot DOPF simulation
figMPC = figure;
plot(t,p_gen_dist,'--d',...
    t,p_demand_dist,'--*',...
    t,bnow_dist,'--o');
hold on;
plot([0 T],[battery_capacity*mpc.Nbus battery_capacity*mpc.Nbus],':k','Linewidth',1);
hold off;
title(strcat('Distributed MPC Simulation'));
xlabel('Time (Hour)');
ylabel('Power (p.u.)');
grid on;
set(gca,'Xtick',0:2:T);
xlim([0 T]);
%ylim([0 0.8]);
legend('Sub Station Generation','Total Demand','Total Storage','Location','southeast');
figstr = 'consensus_sim_MPC';
pngname = strcat(figstr,'.png');
print(figMPC,pngname,'-dpng');
%}

%% Sim All 

% -------------------------------------
%             Prerequisite 
% 1) load gurobi 
% load('mpc_gurobi_storage_0.35.mat');
% 2) run open loop 
% 3) run post analysis 
% 4) run sim all
% -------------------------------------

% open loop computation 
bOpenLoop = sum(bOffline);              % open loop control take the offline battery scheduling
p_b_OpenLoop = bOpenLoop(1,2:end) - bOpenLoop(1,1:end-1);
genOpenLoop = sum(demandOnline*mpc.pf) + [lossOffline lossOffline(1)] + [p_b_OpenLoop p_b_OpenLoop(1)]; 

% plot mpc simulation 
figAll = figure;
t = 0:1:T;
% Offline Demand 
plot(t,-[sum(-demand*mpc.pf) sum(-demand(:,1)*mpc.pf)],'r:o');
% Online Demand 
hold on;
plot(t,sum(demandOnline*mpc.pf),'r-o');         % = p_demand_dist
% Open Loop 
hold on;
plot(t,real(genOpenLoop),'g:x',...
    t,bOpenLoop,'b:x');
% Closed-loop Centralized 
hold on;
plot(t,real(genOnline),'g-s',...
    t,sum(b0Online),'b-s');
% Closed-loop Distributed
hold on;
plot(t,p_gen_dist,'g-p',...
    t,bnow_dist,'b-p');
% Storage Capacity
hold on;
plot([0 T],[battery_capacity*mpc.Nbus battery_capacity*mpc.Nbus],'-.k',...
    [0 T],[0 0],'-.k');
hold off;
legend('Offline Demand','Online Demand',...
    'Open Loop Generation','Open Loop Storage',...
    'MPC Generation Centralized','MPC Storage Centralized',...
    'MPC Generation Centralized','MPC Storage Centralized',...
    'Storage Capacity','Location','east');
title('Demand Response Simulation');
xlabel('Time (Hour)');
ylabel('Power (p.u.)');
set(gca,'Xtick',0:2:T);
xlim([0 T]);
grid on;
set(figAll, 'Position', [0 0 1500 700]);
figstr = 'all_sim_MPC';
pngname = strcat(figstr,'.png');
print(figAll,pngname,'-dpng');


%% Network Balance 

% closed-loop centralized OPF 
powergapCent = real(genOnline) - sum(demandOnline*mpc.pf) - lossOnline; 
bCent(1,1) = sum(b0);
for t=1:T
    sumgapCent(1,t) = sum(powergapCent(1:t));
end
bCent(1,2:T+1) = sum(b0) + sumgapCent(1:T);
unbalCent = sum(b0Online) - bCent;

% closed-loop distributed OPF 
loss_dist = sum(L_dist.*repmat(mpc.R,1,T+1));
powergapDist = p_gen_dist - p_demand_dist - loss_dist'; % p_demand_dist = sum(demandOnline*mpc.pf)'
bDist(1,1) = sum(b0);
for t=1:T
    sumgapDist(1,t) = sum(powergapDist(1:t));
end
bDist(1,2:T+1) = sum(b0) + sumgapDist(1:T);
unbalDist = bnow_dist' - bDist;

% plot unbalance 
figBal = figure;
t = 0:1:T;
plot(t,unbalCent,'--d',...
    t,unbalDist,'--*');
legend('Centralized','Distributed');
grid on;
xlim([0 T]);
set(gca,'Xtick',0:2:T);
ylim([-battery_capacity*mpc.Nbus battery_capacity*mpc.Nbus]);
title('Generation-Storage Balance');
ylabel('Unbalance (p.u.)');
figstr = strcat('bal_MPC');
pngname = strcat(figstr,'.png');
print(figBal,pngname,'-dpng');

%% Warm Starting Analysis 
close all;
% Performance metric  
exactMean = mean(mean(exact_dist));                     
unbalMean = mean(abs(unbalDist));                       
primalStartMean = mean(primalStart);
primalEndMean = mean(primalEnd);
dualStartMean = mean(dualStart);
dualEndMean = mean(dualEnd);
objMisMean = mean(abs(objRatio-1));
batMisMean = mean(abs(bnowRatio-1));
genMisMean = mean(abs(genRatio-1));

disp('---------------------------------');
disp('Status: Well done, oh yeah !!!');
disp('=================================');

cd(rootDir);
resultfile = strcat(TypeAlg,'_rho_',num2str(RHO),'_iter_',num2str(Niter),'_pert_',num2str(PERTMIN),'-',num2str(PERTMAX),'.mat');
if exist(resultfile,'file')
    delete(resultfile)
end
save(resultfile);

set(0,'DefaultFigureVisible','on');

%end