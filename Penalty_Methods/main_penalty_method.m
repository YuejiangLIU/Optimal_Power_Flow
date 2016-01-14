% penalty method
function main()

clear all;  clc;

mu = 0.1;
beta = 10;
epsilon = 1e-3;

x0 = [0;0];

% -----------------------------------
% obj: f(x) = x(1)^2 + 2*x(2)^2;
% const: g(x) = 1 - x1 - x2 == 0;
% penaly function: p(x) = (1-x1-x2)^2;
% h(x) = f(x) + mu * p(x)
% ------------------------------------

x = x0;
res = 1;
niter = 0;
fprintf('iter\tres\t\t\tx1\t\tx2\t\tmu\n');  
while res > epsilon & niter < 10
    % minimization 
    x = [ 2*mu/(2+3*mu) ; mu/(2+3*mu) ];
    % residual
    res = mu * pfun(x);
    % disp 
    niter = niter + 1;
    fprintf('%d\t\t%.4f\t\t%.2f\t%.2f\t%.1f\n',niter,res,x(1),x(2),mu);
    % penalty parameter
    mu = mu * beta;
end

end

function p = pfun(x)
    p = (1-x(1)-x(2))^2;
end