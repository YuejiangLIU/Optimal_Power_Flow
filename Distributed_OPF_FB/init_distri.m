function [xi,zi,di] = init_distri(Nele,Ngen,Tpred)
% X Initialization 
xi = ones(Nele,1);

% Z Initialization 
zi = ones(4,1);

% Dual Initialization 
di.lambda = ones(Tpred*3,1);
di.gamma = ones(Tpred,1);
di.mu = ones(Nele-Ngen*2*Tpred,1);

end