function f = FunF(mat,z,dual,RHO)

Ai = mat.A;
bi = mat.b;
Di = mat.D;
Ei = mat.E;
Fi = mat.F;

lambda = dual.lambda;
gamma = dual.gamma;
mu = dual.mu;

f = @lagrangian_opf;
    function v = lagrangian_opf(x)
        v = Fi'*x + lambda'*(Ai*x-bi) + RHO/2*(Ai*x-bi)'*(Ai*x-bi) ...
            + mu'*(Di*x-z) + RHO/2*(Di*x-z)'*(Di*x-z) ... 
            + gamma'*(x'*Ei*x) + RHO/2*(x'*Ei*x)'*(x'*Ei*x);
    end

end