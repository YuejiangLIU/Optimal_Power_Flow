function xi_new = distr_x_opt(xi,mat,zi_hat,di,RHO)
%:: min Fi'*xi 
%::     + lambda_i'*(Ai*xi - bi) + rho/2*norm(Ai*xi-bi)^2 
%::     + gamma_i'*(xi'*Ei*xi) + rho/2*norm(xi'*Ei*xi)^2 
%::     + mu_i'*(Di*xi - zi_hat) + rho/2*norm(Di*xi - zi_hat)
%:: s.t. LBi <= xi <= UBi

INIT_Ci = 1;
KMAX = 2;
NMAX = 10;
EPSILON = 1e-4;
k = 1;
n = 1;
beta = 2;

% root bus xi(2:4) == 0 
if mat.E(1,2) < 0.1 
    xi(2:4) = 0;
end 

funf = FunF(mat,zi_hat,di,RHO);
proxg = ProxG(mat.LB,mat.UB);
gradf = GradF(mat,zi_hat,di,RHO);

while k <= KMAX
    %{
    Li = RHO * (mat.A'*mat.A + mat.D'*mat.D);
    Li(1:4,1:4) = Li(1:4,1:4) + ...
              RHO * [ xi(2)^2, 2*xi(1)*xi(2)-xi(3)^2-xi(4)^2, -2*xi(2)*xi(3), -2*xi(2)*xi(4); ...
                    2*xi(1)*xi(2)-xi(3)^2-xi(4)^2, xi(1)^2, -2*xi(1)*xi(3), -2*xi(1)*xi(4); ...
                    -2*xi(2)*xi(3), -2*xi(1)*xi(3), 6*xi(3)^2-2*xi(1)*xi(2)+2*xi(4)^2, 4*xi(3)*xi(4); ...  
                    -2*xi(2)*xi(4), -2*xi(1)*xi(4), 4*xi(3)*xi(4), 6*xi(4)^2-2*xi(1)*xi(2)+2*xi(3)^2 ];
    yi = xi - Li \ gradf(xi);
    xi_new = proxg(yi);
    %}
    
    Ci = INIT_Ci;
    while n <= NMAX
        yi = xi - 1/Ci * gradf(xi);
        xi_new = proxg(yi);
        if funf(xi_new) <= funf(xi) + gradf(xi)'*(xi_new-xi) + Ci/2*(xi_new-xi)'*(xi_new-xi)
            break;
        else
            n = n+1;
            Ci = beta * Ci;
        end
    end
    if n>= NMAX
        %disp('backtracking truncated');
        %disp(k);
        fprintf('Ci = %.2f\n',Ci);
        k
        [xi, xi_new]
        pause(1e-3);
    end
    n = 1;

    % root bus xi(2:4) == 0 
    if mat.E(1,2) < 0.1
        xi_new(2:4) = 0;
    end
    
    if max(abs(xi_new./xi-1)) <= EPSILON
        xi = xi_new;
        break;
    else
        k = k+1;
        xi = xi_new;
    end
end

%{
if k >= KMAX 
    disp('x-opt not converge yet');
end
%}

%fprintf('Niter = %d\n',k);
%fprintf('xi = \n',xi);


%% test formulation 
%{
yalmip('clear');
clear xi
% var
xi = sdpvar(length([mat.LB]),1);
% obj 
obj = mat.F'*xi ...
    + di.lambda'*(mat.A*xi-mat.b) + RHO/2*(mat.A*xi-mat.b)'*(mat.A*xi-mat.b) ...
    + di.mu'*(mat.D*xi-zi_hat) + RHO/2*(mat.D*xi-zi_hat)'*(mat.D*xi-zi_hat) ...
    + di.gamma'*(xi'*mat.E*xi) ; % + rho/2*(xi'*mat.E*xi)'*(xi'*mat.E*xi) ...
% constraint
con = [mat.LB <= xi <= mat.UB];
% options
ops = sdpsettings('verbose',0,'solver','gurobi');
% solve
sol = optimize(con,obj,ops);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution = value(xi)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
%}

%{
yalmip('clear');
clear xi
% var
xi = sdpvar(length([mat.LB]),1);
% obj 
obj = mat.F'*xi + di.mu'*(mat.D*xi-zi_hat) + RHO/2*(mat.D*xi-zi_hat)'*(mat.D*xi-zi_hat);
% constraints
con = [mat.LB <= xi <= mat.UB];
con = con + [mat.A * xi == mat.b];
% con = con + [xi' * mat.E * xi >= 0];
if mat.E(1,2) > 0
    con = con + [ norm([2*xi(3);2*xi(4);xi(1)-xi(2)],2) <= xi(1)+xi(2) ];
else
    con = con + [xi(2:4) == 0];
end
% options
ops = sdpsettings('verbose',0,'solver','gurobi');
% solve
sol = optimize(con,obj,ops);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 xi_new = double(xi);
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end
%}

end