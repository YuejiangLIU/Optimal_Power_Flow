%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       autorun a bunch of distopf simulations   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all
clc;

PertList = 5:7;
PertList(PertList == 4) = [];
NiterList = 500;
RhoList = 0.1;

for iRho = 1:length(RhoList)
    RHO = RhoList(iRho);
    for iIter = 1:length(NiterList)
        Niter = NiterList(iIter);
        for iPert = 1:length(PertList)
            PertUnit = PertList(iPert);
            distopt(PertUnit,Niter,RHO);
        end
    end
end
return;
