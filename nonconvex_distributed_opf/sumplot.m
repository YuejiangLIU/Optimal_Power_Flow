function [] = sumplot(off_mpc,net_mpc,sol_opf_mpc,sol_mpc_ref,sol_socp_mpc)

% :: summary plot

tt = 1:net_mpc.Tpred;

figure;
subplot(3,1,1);
plot(tt,sum(net_mpc.p_demand),'-o')
hold on 
plot(tt,sum(off_mpc.p_demand),'--')
ylim([min(min(sum(net_mpc.p_demand),sum(off_mpc.p_demand)))*0.9 max(max(sum(net_mpc.p_demand),sum(off_mpc.p_demand)))*1.1 ]);
legend('online real demand','offline predicted demand','Location','northwest');
ylabel('Demand (p.u.)');
subplot(3,1,2);
plot(tt,sum(sol_opf_mpc.p_gen,1)','-o');
hold on
plot(tt,sum(off_mpc.p_gen,1),'--');
legend('Algorithm 2','offline scheduling','Location','northwest');
ylabel('Generation (p.u.)');


subplot(3,1,3);
plot([tt net_mpc.Tpred+1],[sum(net_mpc.e_bat_now,1) sum(sol_opf_mpc.e_bat,1)],'r-o');
hold on
plot([tt net_mpc.Tpred+1],[sum(net_mpc.e_bat_now,1) sum(sol_mpc_ref.e_bat,1)],'b-s');
hold on
plot([tt net_mpc.Tpred+1],[sum(net_mpc.e_bat_now) sum(off_mpc.e_bat,1)],'g--*');
ylabel('Battery Energy (p.u.)');
legend('Algorithm 2','closed-loop optimal','offline scheduling');


figure;
subplot(2,1,1);
plot(tt,(sol_opf_mpc.obj./sol_mpc_ref.obj - 1),'o-');
hold on;
plot(tt,(sol_socp_mpc.obj./sol_mpc_ref.obj - 1),'k--');
ylabel('Objective Diff');
% ylim([-1e-3 1e-3]);

subplot(2,1,2);
plot(tt,sol_opf_mpc.res,'o-');
hold on;
plot(tt,sol_socp_mpc.res,'k--');
xlabel('Hour');
ylabel('Feasiblity Error');
set(gca,'yscale','log');
% ylim([min(sol_opf_mpc.res)*0.9 max(sol_mpc_ref.res)*1.1]);

end