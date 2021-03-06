clear all;
clc;
close all;

psi = 2;

a = 2;
b = 1;
c1 = 4;
c2 = 2;

Ts = 0.1:0.001:3;
dW = ( c1 * (a * Ts).^(-psi) + c2 ) .* ( b * Ts); 

%%
close all
figure;
plot(Ts,dW);
hold on;
hlb = plot([1/a 1/a],[min(dW)*0.9 max(dW)*1.1],'k--','Linewidth',1.5);
xlabel('Ts');
ylabel('f');
legend(hlb,'Ts lower bound');
% ylim([0 10]);

annotation('textbox',...
    [0.7 0.15 0.15 0.27],...
    'String',{['a = ' num2str(value(a))],...
    ['b = ' num2str(value(b))],...
    ['c1 = ' num2str(value(c1))],...
    ['c2 = ' num2str(value(c2))],...
    ['\Psi = ' num2str(value(psi))]},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','--',...
    'EdgeColor',[1 1 0],...
    'LineWidth',1,...
    'BackgroundColor',[0.9  0.9 0.9],...
    'Color',[0.84 0.16 0]);

[~,idx] = min(dW);
Ts_extreme_est = 1/a * (c1 * (psi-1) / c2)^(1/psi) 
Ts_extreme_sim = (idx-1) * 0.001 + 0.1 