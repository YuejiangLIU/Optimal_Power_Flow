function loadpred = pred(load,npred,plotind,pertmin,pertmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Real-time Distributed Optimal Power Flow 
%       LA @ EPFL, Jul. 26, 2015 
%       Contact: yuejiang.liu@epfl.ch / liuyuejiang1989@gmail.com
%       Prerequisite: YALMIP, Gurobi, (CVX)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');     % generate same perturbation seed 
DOPLOT = 1;
T = size(load,2);
loadshift = load;
for ii = 1:npred
    
    loadpred{ii} = pertfunction(loadshift,pertmin,pertmax);

    % >>> Plot Section Start <<< 
    if DOPLOT && any(ii == plotind)
    figPred = figure;
    stairs([loadpred{ii}';loadpred{ii}(:,end)'],'LineWidth',1);
    hold on;
    stairs([loadshift';loadshift(:,end)'],':','LineWidth',1);
    hold off;
    xlim([1 T+1]);
    xtick = 1:2:T+1;
    set(gca,'XTick',xtick);
    set(gca,'XTickLabel',xtick-2+ii);
    title(['Hourly Demand Prediction at Time ' num2str(ii-1)]);
    xlabel('Time (Hour)');
    ylabel('Load (p.u.)');
    grid on;
    legend('Bus 1','Bus 2','Bus 3','Bus 4','Bus 5','Bus 6','Bus 7','Location','northeast');

    figstr = strcat('pred_time_',num2str(ii-1));
    %figname = strcat(figstr,'.fig');
    %savefig(figname);
    pngname = strcat(figstr,'.png');
    saveas(figPred,pngname,'png');
    close all;
    end
    % >>> Plot Section End <<< 

    loadshift = circshift(loadpred{ii},-1,2);

end

end

function loadpert = pertfunction(load,pertmin,pertmax)            
pertrange = pertmax - pertmin;
Nbus = size(load,1);     
pertabs = rand(Nbus,1)*pertrange + pertmin;
if rand>0.5
    pertdir = 1;            % all buses are perturbed to same direction
else 
    pertdir = -1;
end
pertreal = pertabs*pertdir;
loadpert = load;
loadpert(:,1) = load(:,1).*(1+pertreal);
end