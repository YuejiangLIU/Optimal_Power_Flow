function [] = rtplot(off_mpc,net_mpc,sol_opf_mpc,sol_socp_mpc)

% :: plot real-time closed loop results

tt = 1:net_mpc.Tpred;

h = figure;
subplot(4,1,1);
plot(tt,sum(net_mpc.p_demand),'-o')
hold on 
plot(tt,sum(off_mpc.p_demand),'--')
ylim([min(min(sum(net_mpc.p_demand),sum(off_mpc.p_demand)))*0.9 max(max(sum(net_mpc.p_demand),sum(off_mpc.p_demand)))*1.1 ]);
legend('online real demand','offline predicted demand','Location','southeast');
ylabel('Demand (p.u.)');
xlim([1 net_mpc.Tpred+1]);
set(gca,'XTick',1:2:net_mpc.Tpred+1);
set(gca,'XTickLabel',sprintf('%d\n',0:2:net_mpc.Tpred));
ylim([3.5 5.8]);

subplot(4,1,2);
plot(tt,sum(sol_opf_mpc.p_gen,1)','-o');
hold on
plot(tt,sum(off_mpc.p_gen,1),'--');
legend('Algorithm 2','offline scheduling','Location','southeast');
ylabel('Generation (p.u.)');
xlim([1 net_mpc.Tpred+1]);
set(gca,'XTick',1:2:net_mpc.Tpred+1);
set(gca,'XTickLabel',sprintf('%d\n',0:2:net_mpc.Tpred));
ylim([3.5 5.8]);

subplot(4,1,3);
plot([tt net_mpc.Tpred+1],[sum(net_mpc.e_bat_now,1) sum(sol_opf_mpc.e_bat,1)],'-o');
hold on 
plot([tt net_mpc.Tpred+1],[sum(net_mpc.e_bat_now,1) sum(off_mpc.e_bat,1)],'--');
ylabel('Battery Energy (p.u.)'); 
legend('Algorithm 2','offline scheduling');
xlim([1 net_mpc.Tpred+1]);
set(gca,'XTick',1:2:net_mpc.Tpred+1); 
set(gca,'XTickLabel',sprintf('%d\n',0:2:net_mpc.Tpred));
ylim([0 0.6]);

subplot(4,1,4);
plot(tt,sol_opf_mpc.res,'o-');
hold on;
plot(tt,sol_socp_mpc.res,'^-');
xlabel('Hour');
ylabel('Feasiblity Error');
set(gca,'yscale','log');
ylim([min(sol_opf_mpc.res)*0.9 max(sol_socp_mpc.res)*4]);
legend('Algorithm 2','convex relaxation');
xlim([1 net_mpc.Tpred+1]);
set(gca,'XTick',1:2:net_mpc.Tpred+1);
set(gca,'XTickLabel',sprintf('%d\n',0:2:net_mpc.Tpred));
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1 10 1e2]);

pos = get(h,'Position');
pos(1:2) = [0 0];
pos(4) = 800;
set(h,'Position',pos);
save2pdf('open_vs_close.pdf',h,600);

end