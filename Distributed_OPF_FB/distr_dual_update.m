function di_new = distr_dual_update(xi,mat,zi_hat,di,RHO)

di_new.lambda = di.lambda + RHO * (mat.A*xi - mat.b);
di_new.gamma = di.gamma + RHO * (xi'*mat.E*xi);
di_new.mu = di.mu + RHO * (mat.D*xi - zi_hat);

end