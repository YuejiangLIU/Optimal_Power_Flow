function funf = FunF(idx,net,RESCALE)

% :: compute f_i(x_i)

Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);
Ngen = net.gen(idx,1);
Nbat = net.bat(idx,1);

funf = @lagrangian_opf;
    function f = lagrangian_opf(xi,z,dual,par,RHO)
        f = 0;
        % obj 
        if Ngen > 0
            if net.ncost == 3
                % fcost = sum( net.c_gen(idx,1) * (xi.bus(1+Nabus*3+1,:) * net.baseMVA).^2 + net.c_gen(idx,2) * (xi.bus(1+Nabus*3+1,:) * net.baseMVA) );
                fcost = sum( net.c_gen(idx,1) * (xi.bus(1+Nabus*3+1,:)).^2 );
            elseif net.ncost == 2
                fcost = sum( net.c_gen(idx,1) * (xi.bus(1+Nabus*3+1,:) * net.baseMVA) ); 
            end
            f = f + fcost * RESCALE;
        end
        % con1&2
        temp1 = - par.p_demand(idx,:);
        temp2 = - par.q_demand(idx,:);
        if Nabus > 0 
            temp1 = temp1 - xi.bus(3,:);
            temp2 = temp2 - xi.bus(4,:);
        end
        if Ngen > 0 
            temp1 = temp1 + xi.bus(1+Nabus*3+1,:);
            temp2 = temp2 + xi.bus(1+Nabus*3+2,:);
        end
        if Nbat > 0 
            temp1 = temp1 + xi.bus(1+Nabus*3+Ngen*2+1,:);
        end
        if Ncbus > 0 
            for kC = 1:Ncbus
                idxC = net.cbus(idx,kC+1);
                temp1 = temp1 + xi.cbus(2+(kC-1)*3,:) - xi.cbus(1+(kC-1)*3,:)*net.ri(idxC);
                temp2 = temp2 + xi.cbus(3+(kC-1)*3,:) - xi.cbus(1+(kC-1)*3,:)*net.xi(idxC);
            end
        end
        f = f + dual{idx}.lambda(1,:) * temp1' + RHO/2 * norm(temp1,2)^2;
        f = f + dual{idx}.lambda(2,:) * temp2' + RHO/2 * norm(temp2,2)^2;
        % con3 
        if Nabus > 0
             temp3 = xi.bus(1,:) + xi.bus(2,:)*(net.ri(idx)^2+net.xi(idx)^2) ...
                - 2*(xi.bus(3,:)*net.ri(idx)+xi.bus(4,:)*net.xi(idx)) - xi.abus; 
             f = f + dual{idx}.lambda(3,:) * temp3' + RHO/2 * norm(temp3,2)^2; 
        end        
        % con battery 
        if Nbat > 0
            tempbat(1,1) = net.e_bat_now(idx) - xi.bus(1+Nabus*3+Ngen*2+1,1) - xi.bus(1+Nabus*3+Ngen*2+2,1);    % e0 - p1 = e1
            if net.Tpred > 1
                for tt = 2:net.Tpred
                    tempbat(1,tt) = xi.bus(1+Nabus*3+Ngen*2+2,tt-1) - xi.bus(1+Nabus*3+Ngen*2+1,tt) - xi.bus(1+Nabus*3+Ngen*2+2,tt);
                end
                f = f + dual{idx}.lambda(2+Nabus+1,:) * tempbat' + RHO/2 * norm(tempbat,2)^2;
            end
        end
        % non-convex con
        if Nabus > 0
            tempnc = xi.bus(1,:) .* xi.bus(2,:) - xi.bus(3,:).^2 - xi.bus(4,:).^2;
            f = f + dual{idx}.gamma * tempnc' + RHO/2 * norm(tempnc,2)^2;
        end
        % consensus con 
        f = f + sum(sum( dual{idx}.mu.bus .* ( xi.bus(1:1+Nabus*3,:) - z{idx}(1:1+Nabus*3,:) ) ))...
            + RHO/2 * norm(xi.bus(1:1+Nabus*3,:) - z{idx}(1:1+Nabus*3,:),2)^2;
        if Nabus > 0 
            idxA = net.abus(idx,2);
            f = f + dual{idx}.mu.abus * (xi.abus - z{idxA}(1,:))' ...
                + RHO/2 * norm(xi.abus - z{idxA}(1,:),2)^2;
        end 
        if Ncbus > 0
            for kC = 1:Ncbus
                idxC = net.cbus(idx,kC+1);
                ztemp(1+(kC-1)*3:3+(kC-1)*3,:) = z{idxC}(2:4,:);
            end
            f = f + sum(sum( dual{idx}.mu.cbus .* (xi.cbus - ztemp) )) ...
                + RHO/2 * norm(xi.cbus - ztemp,2)^2;
        end
    end
end