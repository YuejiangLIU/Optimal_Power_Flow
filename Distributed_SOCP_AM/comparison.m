%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       result analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
%% comparison of perturbation 
clear all;
close all
clc;

PertList = 1:4;
NiterList = 500;
RhoList = 0.1;

for iRho = 1:length(RhoList)
    RHO = RhoList(iRho);
    for iIter = 1:length(NiterList)
        Niter = NiterList(iIter);
        for iPert = 1:length(PertList)
            PertUnit = PertList(iPert);
            matname = strcat('ADMM_warm_rho_',num2str(RHO),'_iter_',num2str(Niter),'_pert_',num2str(PertUnit),'-',num2str(PertUnit),'.mat');
            load(matname,'exact_dist','unbalADMM','powergapADMM','bnow_dist','primalStartMean','primalEndMean','dualStartMean','dualEndMean','objMisMean','batMisMean','genMisMean');
            %exactMaxList(iRho) = exactMax;
            exactMeanList(iPert,1) = mean(mean(exact_dist(2:end,:)));
            %unbalMeanList(iRho,iIter) = mean(abs(powergapADMM));
            p_battery_dist = bnow_dist(2:end) - bnow_dist(1:end-1);
            unbalMeanList(iPert,1) = mean(abs(powergapADMM(1:24) - p_battery_dist));
            %unbalMeanList(iPert) = unbalADMM(end);
            primalStartMeanList(iPert,1) = primalStartMean;
            primalEndMeanList(iPert,1) = primalEndMean;
            dualStartMeanList(iPert,1) = dualStartMean;
            dualEndMeanList(iPert,1) = dualEndMean;
            objMisList(iPert,1) = objMisMean;
            batMisList(iPert,1) = batMisMean;
            genMisList(iPert,1) = genMisMean;
        end
    end
end

figVal = figure;
subplot(3,1,1);
plot(PertList,exactMeanList(:,1),'-o');
grid on;
ylim([1 1.01]);
ylabel('Exact Mean');
title('ADMM Solution');
subplot(3,1,2);
plot(PertList,unbalMeanList(:,1),'-o');
grid on;
ylabel('Unbalance Mean (p.u.)');
subplot(3,1,3);
plot(PertList,objMisList(:,1),'-o');
ylabel('Objective Diff');
xlabel('Perturbation (%)');
%title('Average Difference ADMM vs Gurobi (%)');
ylim([0 0.01]);
grid on;
set(figVal, 'Position', [0 0 1500 800]);
%legend('ADMM Niter = 100','Homotopy = 10 Niter = 500');
print(figVal,'ADMM_sol_cmp_pert.png','-dpng');

figRes = figure;
plot(PertList,primalStartMeanList(:,1),'-o',PertList,primalEndMeanList(:,1),'-o',...
    PertList,dualStartMeanList(:,1),'-.o',PertList,dualEndMeanList(:,1),'-.o');
set(gca,'yscale','log');
title('ADMM Residual');
grid on;
legend('Primal Residual Start','Primal Residual End','Dual Residual Start','Dual Residual End','Location','northeast');
xlabel('Perturbation (%)');
set(figRes, 'Position', [0 0 1200 700]);
print(figRes,'ADMM_res_cmp_pert.png','-dpng');

return;
figObj = figure;
subplot(2,1,1);
plot(PertList,batMisList,'-o');
ylabel('Storage');
ylim([0 0.02]);
title('ADMM Input');
grid on;
subplot(2,1,2);
plot(PertList,genMisList,'-o');
ylabel('Generator');
ylim([0 0.02]);
grid on;
xlabel('Perturbation');
set(figObj, 'Position', [0 0 1200 700]);
print(figObj,'ADMM_input_cmp_pert.png','-dpng');


%% comparison of iteration 
clear all;
close all
clc;

PertList = 4;
NiterList = 100:100:700;
RhoList = 0.1;

for iRho = 1:length(RhoList)
    RHO = RhoList(iRho);
    for iIter = 1:length(NiterList)
        Niter = NiterList(iIter);
        for iPert = 1:length(PertList)
            PertUnit = PertList(iPert);
            matname = strcat('ADMM_warm_rho_',num2str(RHO),'_iter_',num2str(Niter),'_pert_',num2str(PertUnit),'-',num2str(PertUnit),'.mat');
            load(matname,'exact_dist','unbalADMM','primalStartMean','primalEndMean','dualStartMean','dualEndMean','objMisMean','batMisMean','genMisMean');
            exactMeanList(iIter) = mean(mean(exact_dist(2:end,:)));
            unbalMeanList(iIter) = mean(abs(unbalADMM));
            primalStartMeanList(iIter) = primalStartMean;
            primalEndMeanList(iIter) = primalEndMean;
            dualStartMeanList(iIter) = dualStartMean/RHO;
            dualEndMeanList(iIter) = dualEndMean/RHO;
            objMisList(iIter) = objMisMean;
            batMisList(iIter) = batMisMean;
            genMisList(iIter) = genMisMean;
        end
    end
end

figVal = figure;
subplot(2,1,1);
plot(NiterList,exactMeanList,'-o');
grid on;
ylim([1 1.01]);
ylabel('Exact Mean');
title('Solution Validation');
subplot(2,1,2);
plot(NiterList,unbalMeanList,'-o');
grid on;
ylabel('Power Unbalance Mean (p.u.)');
xlabel('Iteratoin Limit');
print(figVal,'ADMM_valid_cmp_niter.png','-dpng');

figRes = figure;
plot(NiterList,primalStartMeanList,'-o',NiterList,primalEndMeanList,'-o',...
    NiterList,dualStartMeanList,'-.o',NiterList,dualEndMeanList,'-.o');
set(gca,'yscale','log');
title('Residual');
grid on;
legend('Primal Residual Start','Primal Residual End','Dual Residual Start','Dual Residual End','Location','northeast');
xlabel('Iteratoin Limit');
print(figRes,'ADMM_res_cmp_niter.png','-dpng');

figObj = figure;
subplot(3,1,1);
plot(NiterList,objMisList,'-o');
ylabel('Objective');
title('Average Difference ADMM vs Gurobi (%)');
ylim([0 0.01]);
grid on;
subplot(3,1,2);
plot(NiterList,batMisList,'-o');
ylabel('Storage');
ylim([0 0.02]);
grid on;
subplot(3,1,3);
plot(NiterList,genMisList,'-o');
ylabel('Generator');
ylim([0 0.02]);
grid on;
xlabel('Iteratoin Limit');
print(figObj,'ADMM_obj_cmp_niter.png','-dpng');


%% comparison of penalty parameter  
clear all;
close all
clc;

PertList = 2;
NiterList = [500];
%RhoList = [1e-3 2e-3 5e-3 1e-2 2e-2 5e-2 1e-1 2e-1 5e-1 1 2 5 10 20 50 100 2e2 5e2 1e3];
%RhoList = [0.0001 0.0005 0.001 0.002 0.004 0.005 0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.5 0.75 1 2 5 10 20 50 60 70 80 100 120 150 200 300 500];
%RhoList = [0.01 0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.5];
RhoList = [1e-2 1e-1 1];
for iRho = 1:length(RhoList)
    RHO = RhoList(iRho);
    for iIter = 1:length(NiterList)
        Niter = NiterList(iIter);
        for iPert = 1:length(PertList)
            PertUnit = PertList(iPert);
            matname = strcat('ADMM_warm_rho_',num2str(RHO),'_iter_',num2str(Niter),'_pert_',num2str(PertUnit),'-',num2str(PertUnit),'.mat');
            load(matname,'exact_dist','powergapADMM','bnow_dist','primalStartMean','primalEndMean','dualStartMean','dualEndMean','objMisMean','batMisMean','genMisMean');
            %exactMaxList(iRho) = exactMax;
            exactMeanList(iRho,iIter) = mean(mean(exact_dist(2:end,:)));
            %unbalMeanList(iRho,iIter) = mean(abs(powergapADMM));
            p_battery_dist = bnow_dist(2:end) - bnow_dist(1:end-1);
            unbalMeanList(iRho,iIter) = mean(abs(powergapADMM(1:24) - p_battery_dist));
            primalStartMeanList(iRho,iIter) = primalStartMean;
            primalEndMeanList(iRho,iIter) = primalEndMean;
            dualStartMeanList(iRho,iIter) = dualStartMean;
            dualEndMeanList(iRho,iIter) = dualEndMean;
            objMisList(iRho,iIter) = objMisMean;
            batMisList(iRho,iIter) = batMisMean;
            genMisList(iRho,iIter) = genMisMean;
        end
    end
end

figVal = figure;
subplot(3,1,1);
plot(RhoList,exactMeanList(:,1),'r-o');
hold off;
grid on;
set(gca,'xscale','log');
xlim([1e-2 1]);
ylim([1 1.01]);
ylabel('Exact Mean');
legend('ADMM Niter = 500','Location','northwest');
title('ADMM Solution');
subplot(3,1,2);
plot(RhoList,unbalMeanList(:,1),'r-o');
hold off;
grid on;
set(gca,'xscale','log');
ylabel('Unbalance Mean (p.u.)');
xlim([1e-2 1]);
ylim([0 2e-2]);
xlabel('Rho');
%legend('Truncated Niter = 100','Truncated Niter = 200','Truncated Niter = 300','Truncated Niter = 500');
subplot(3,1,3);
plot(RhoList,objMisList(:,1),'r-o');
hold off;
ylabel('Objective Diff');
%title('Average Difference ADMM vs Gurobi (%)');
xlim([1e-2 1]);
ylim([0 0.01]);
grid on;
set(gca,'xscale','log');
%legend('Truncated Niter = 100','Truncated Niter = 200','Truncated Niter = 300','Truncated Niter = 500');
set(figVal, 'Position', [0 0 1500 800]);
figValName = strcat('ADMM_sol_cmp_rho_niter.png');
print(figVal,figValName,'-dpng');

figRes = figure;
plot(RhoList,primalStartMeanList(:,1),'r-^',RhoList,primalEndMeanList(:,1),'r-v',...
    RhoList,dualStartMeanList(:,1),'r-.o',RhoList,dualEndMeanList(:,1),'r-.*');
hold off;
set(gca,'yscale','log');
title('ADMM Residual');
grid on;
set(gca,'xscale','log');
legend('Primal Residual Start','Primal Residual End','Dual Residual Start','Dual Residual End','Location','northeast');
xlim([1e-2 1]);
xlabel('Rho');
%legend('Niter = 100','Niter = 200','Niter = 300');
%figResName = strcat('ADMM_res_cmp_rho_niter',num2str(Niter),'.png');
set(figRes, 'Position', [0 0 1200 700]);
figResName = strcat('ADMM_res_cmp_rho_niter.png');
print(figRes,figResName,'-dpng');

return;

figObj = figure;
subplot(2,1,1);
plot(RhoList,batMisList(:,1),'r-o');
hold on;
plot(RhoList,batMisList(:,2),'g-o');
hold on;
plot(RhoList,batMisList(:,3),'b-o');
hold on;
plot(RhoList,batMisList(:,4),'k-o');
hold off;
ylabel('Storage');
xlim([1e-3 1e3]);
ylim([0 0.02]);
title('ADMM Input');
grid on;
%legend('Truncated Niter = 100','Truncated Niter = 200','Truncated Niter = 300','Truncated Niter = 500');
set(gca,'xscale','log');
subplot(2,1,2);
plot(RhoList,genMisList(:,1),'r-o');
hold on;
plot(RhoList,genMisList(:,2),'g-o');
hold on;
plot(RhoList,genMisList(:,3),'b-o');
hold on;
plot(RhoList,genMisList(:,4),'k-o');
hold on;
ylabel('Generator');
xlim([1e-3 1e3]);
ylim([0 0.02]);
grid on;
set(gca,'xscale','log');
xlabel('Rho');
legend('ADMM Niter = 100','ADMM Niter = 200','ADMM Niter = 300','ADMM Niter = 500','Location','northwest');
set(figObj, 'Position', [0 0 1200 700]);
figObjName = strcat('ADMM_input_cmp_rho_niter.png');
print(figObj,figObjName,'-dpng');

