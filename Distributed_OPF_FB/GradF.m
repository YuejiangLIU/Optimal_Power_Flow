function df = GradF(mat,z,dual,RHO)

Ai = mat.A;
bi = mat.b;
Di = mat.D;
Ei = mat.E;
Fi = mat.F;

lambda = dual.lambda;
gamma = dual.gamma;
mu = dual.mu;

df = @gradient_opf;
    function v = gradient_opf(x)
        v = Fi + Ai'*lambda + RHO*Ai'*(Ai*x-bi) ...
            + Di'*mu + RHO*Di'*(Di*x-z);
        v(1:4,1) = v(1:4,1) + gamma * [x(2,1);x(1,1);-2*x(3,1);-2*x(4,1)] ...
                  + RHO * [x(1,1)*x(2,1)^2 - x(2,1)*x(3,1)^2 - x(2,1) * x(4,1)^2; ...
                          x(1,1)^2*x(2,1) - x(1,1)*x(3,1)^2 - x(1,1)*x(4,1)^2; ...
                          2*x(3,1)^3 - 2*x(1,1)*x(2,1)*x(3,1) + 2*x(3,1)*x(4,1)^2; ... 
                          2*x(4,1)^3 - 2*x(1,1)*x(2,1)*x(4,1) + 2*x(3,1)^2*x(4,1)];  
    end

end