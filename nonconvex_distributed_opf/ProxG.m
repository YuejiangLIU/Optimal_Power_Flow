function proxg = ProxG(idx,net)

% % :: compute projection to omega_i, modified from cvxgrp@github 

Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);
Ngen = net.gen(idx,1);
Nbat = net.bat(idx,1);

proxg = @projection_box; 
    function pi = projection_box(xi,t0) 
        pi.bus(1,:) = max(net.v_min(idx)^2, min(xi.bus(1,:), net.v_max(idx)^2));
        if Nabus > 0 
            pi.bus(2:4,:) = xi.bus(2:4,:);
        end
        if Ngen > 0
            pi.bus(1+Nabus*3+1,:) = max(net.p_gen_min(idx), min(xi.bus(1+Nabus*3+1,:), net.p_gen_max(idx)));
            % pi.bus(1+Nabus*3+1,:) = net.p_gen_ext(idx,:);
            pi.bus(1+Nabus*3+2,:) = max(net.q_gen_min(idx), min(xi.bus(1+Nabus*3+2,:), net.q_gen_max(idx)));
        end
        if Nbat > 0 
            pi.bus(1+Nabus*3+Ngen*2+1,:) = max(net.p_bat_min(idx), min(xi.bus(1+Nabus*3+Ngen*2+1,:), net.p_bat_max(idx)));
            pi.bus(1+Nabus*3+Ngen*2+2,:) = max(net.e_bat_min(idx), min(xi.bus(1+Nabus*3+Ngen*2+2,:), net.e_bat_max(idx)));
            pi.bus(1+Nabus*3+Ngen*2+2,net.Tpred) = net.e_bat_now(idx);
            if net.Tpred == 24
                pi.bus(1+Nabus*3+Ngen*2+2,24-t0) = net.e_bat_midnight(idx);
            end
        end
        if Nabus > 0 
            idxA = net.abus(idx,2);
            pi.abus = max(net.v_min(idxA)^2, min(xi.abus, net.v_max(idxA)^2));
        end
        if Ncbus > 0 
            pi.cbus = xi.cbus;
        end
    end
end