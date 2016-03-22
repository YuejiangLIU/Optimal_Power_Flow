function [] = iterplot(sol,sol_ref)

% :: plot res & obj against iterations 

Niter = length(sol.obj);

res = sqrt(sol.resPF.^2 + sol.resPrim.^2);
obj = sol.obj; 

res_ref = norm(sol_ref.v .*sol_ref.L - (sol_ref.P.^2 + sol_ref.Q.^2));

h = figure;
xx = 1:Niter;
subplot(2,1,1);
hRes = plot(xx,res);
ylim([1e-4 1]);
ylabel('Feasibility Error');
set(gca,'yscale','log');
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1]);

subplot(2,1,2);
hObj = plot(xx,obj/sol_ref.obj-1);
hold on;
hObjRef = plot([1,xx(end)],[0,0],'k--');
ylim([-0.04 0.12]);
set(gca,'YTick',[-0.04 0 0.04 0.08 0.12]);
ylabel('Relative Objective Diff');
xlabel('Iteration Number');
drawnow;

pos = get(h,'Position');
pos(1:2) = [0 0];
pos(4) = 400;
set(h,'Position',pos);
save2pdf('9bus_24h_conv.pdf',h,600);

end