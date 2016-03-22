function [p_pert,q_pert] = loadpred(net,pertmax)

% :: generate online demand profiles by perturbing the receding horizon
% demand prediction  

DOPLOT = 0;

if net.Ncalc >= 1 
    for ii = 1:net.Ncalc  
        p_pert{ii} = pertfun(net.p_demand,ii,pertmax);
        q_pert{ii} = pertfun(net.q_demand,ii,pertmax);
    end
end

end

function d_pert = pertfun(d_base,idx,pertmax)            

d_pert = circshift(d_base,-idx+1,2);

Nbus = size(d_base,1);
pertmag = rand*pertmax;
if rand > 0.7 				% tune random direction probability  
    pertdir = 1;            % all buses are perturbed to same direction
else
    pertdir = -1;
end
pertreal = pertmag*pertdir;

d_pert(:,1) = d_pert(:,1)*(1+pertreal);

end