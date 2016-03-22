function xi_new = distr_x_opt(idx,xi,z,dual,net,par,RHO,t0,funf,gradf,proxg)

% :: x_i update 


INIT_Ci = 1;
KMAX = 1;
NMAX = 50;
EPSILON = 1e-4;

k = 1;

beta = 2;
theta = 0;

Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);
Ngen = net.gen(idx,1);
Nbat = net.bat(idx,1);

while k <= KMAX
    Ci = INIT_Ci;
    n = 1;
    while n <= NMAX

        dxi = gradf(xi,z,dual,par,RHO);
        
        yi.bus = xi.bus - 1/Ci * dxi.bus;
        if Nabus > 0
            yi.abus = xi.abus - 1/Ci * dxi.abus;
        end
        if Ncbus > 0
            yi.cbus = xi.cbus - 1/Ci * dxi.cbus;
        end
      
        xi_new = proxg(yi,t0);
        
        temptaylor = sum(sum( dxi.bus.*(xi_new.bus - xi.bus) + Ci/2*(xi_new.bus-xi.bus).*(xi_new.bus-xi.bus) ));
        if Nabus > 0 
            temptaylor = temptaylor + dxi.abus*(xi_new.abus - xi.abus)' + Ci/2*(xi_new.abus-xi.abus)*(xi_new.abus-xi.abus)';
        end
        if Ncbus > 0
            temptaylor = temptaylor + sum(sum( dxi.cbus.*(xi_new.cbus - xi.cbus) + Ci/2*(xi_new.cbus-xi.cbus).*(xi_new.cbus-xi.cbus) ));
        end  
   
        if funf(xi_new,z,dual,par,RHO) < funf(xi,z,dual,par,RHO) + temptaylor
            break;
        else
            n = n+1;
            Ci = beta * Ci;
        end
       
    end
    
    if n >= NMAX 
        % disp('backtracking not converge yet');
        xi_new = xi;
    end
    
    k = k+1;
    xi = xi_new;

end

end