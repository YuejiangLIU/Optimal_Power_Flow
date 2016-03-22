function simplot(net,sol,t0)

% :: plot predicted solution at time step t0

t = 1:net.Tpred;

plot(t,sum(sol.p_gen,1),'--d',...
    t,sum(net.p_demand,1),'--*',...
    [t,t(end)+1],[sum(net.e_bat_now,1) sum(sol.e_bat,1)],'--o');
xlabel('Time (Hour)');
ylabel('Power (p.u.)');
set(gca,'XTick',1:2:net.Tpred+1);
set(gca,'XTickLabel',sprintf('%d\n',t0:2:net.Tpred+t0));
grid on;
legend('Generation','Demand','Storage','Location','southeast');

p_loss = net.ri'*sol.L;

p_res = sum(sol.p_gen,1) + sum(sol.p_bat,1) - sum(net.p_demand,1) - sum(p_loss,1);

fprintf('Power balance = %.4f\n',sum(p_res));

end