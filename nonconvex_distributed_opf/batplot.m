function [] = batplot(sol,sol_ref)

figure;
xx = 1:length(sol.e_bat);
plot(xx,sol.e_bat);
hold on;
plot(xx,sum(sol_ref.e_bat,1));
legend('sol','sol_ref');

end