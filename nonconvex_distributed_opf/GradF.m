function gradf = GradF(idx,net,RESCALE)

% :: compute df_i(x_i)

Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);
Ngen = net.gen(idx,1);
Nbat = net.bat(idx,1);

gradf = @gradient_opf;
    function df = gradient_opf(xi,z,dual,par,RHO)
        temp1 = - par.p_demand(idx,:);
        temp2 = - par.q_demand(idx,:);
        if Nabus > 0 
            temp1 = temp1 - xi.bus(3,:);
            temp2 = temp2 - xi.bus(4,:);
            temp3 = xi.bus(1,:)+(net.ri(idx)^2+net.xi(idx)^2)*xi.bus(2,:) - 2*(net.ri(idx)*xi.bus(3,:)+net.xi(idx)*xi.bus(4,:))-xi.abus;
            tempnc = xi.bus(1,:).*xi.bus(2,:)-xi.bus(3,:).^2-xi.bus(4,:).^2;
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
                temp1 = temp1 + xi.cbus(2+(kC-1)*3,:) - xi.cbus(1+(kC-1)*3,:)*net.ri(idxC,1);
                temp2 = temp2 + xi.cbus(3+(kC-1)*3,:) - xi.cbus(1+(kC-1)*3,:)*net.xi(idxC,1);
            end
        end
        % local bus gradient 
        df.bus(1,:) = dual{idx}.mu.bus(1,:) + RHO * ( xi.bus(1,:)-z{idx}(1,:) ) ;
        if Nabus > 0 
            df.bus(1,:) = df.bus(1,:) + dual{idx}.lambda(3,:) + RHO * temp3 ...
                + dual{idx}.gamma .* xi.bus(2,:) + RHO * tempnc .* xi.bus(2,:);
        end
        
        if Nabus > 0
            df.bus(2,:) = dual{idx}.mu.bus(2,:) + RHO * ( xi.bus(2,:)-z{idx}(2,:) ) ...
                + dual{idx}.lambda(3,:) * (net.ri(idx)^2+net.xi(idx)^2) + RHO * temp3 * (net.ri(idx)^2+net.xi(idx)^2) ... 
                + dual{idx}.gamma .* xi.bus(1,:) + RHO * tempnc .* xi.bus(1,:);
        
            df.bus(3,:) = dual{idx}.mu.bus(3,:) + RHO * ( xi.bus(3,:)-z{idx}(3,:) ) ...
                - 2 * net.ri(idx) * dual{idx}.lambda(3,:) + RHO * temp3 * (-2*net.ri(idx)) ...
                + dual{idx}.gamma .* (-2*xi.bus(3,:)) + RHO * tempnc .* (-2*xi.bus(3,:)) ...
                - dual{idx}.lambda(1,:) - RHO * temp1;
            
            df.bus(4,:) = dual{idx}.mu.bus(4,:) + RHO * (xi.bus(4,:)-z{idx}(4,:)) ... 
                - 2 * net.xi(idx) * dual{idx}.lambda(3,:) + RHO * temp3 * (-2*net.xi(idx)) ... 
                + dual{idx}.gamma .* ( -2*xi.bus(4,:) ) + RHO * tempnc .* (-2*xi.bus(4,:)) ... 
                - dual{idx}.lambda(2,:) - RHO * temp2;
        end
        
        if Ngen > 0 
            if net.ncost == 3
                % df.bus(1+Nabus*3+1,:) = ( 2 * net.c_gen(idx,1) * xi.bus(1+3*Nabus+1,:) * net.baseMVA^2 + net.c_gen(idx,2) * net.baseMVA ) * RESCALE + dual{idx}.lambda(1,:) + RHO * temp1;
                df.bus(1+Nabus*3+1,:) = ( 2 * net.c_gen(idx,1) * xi.bus(1+3*Nabus+1,:) ) * RESCALE + dual{idx}.lambda(1,:) + RHO * temp1;
            elseif net.ncost == 2
                df.bus(1+Nabus*3+1,:) = net.c_gen(idx,1) * net.baseMVA * RESCALE + dual{idx}.lambda(1,:) + RHO * temp1;
            end
            df.bus(1+Nabus*3+2,:) = dual{idx}.lambda(2,:) + RHO * temp2;               
        end
        
        if Nbat > 0 
            tempbat(1,1) = net.e_bat_now(idx) - xi.bus(1+Nabus*3+Ngen*2+1,1) - xi.bus(1+Nabus*3+Ngen*2+2,1);    % e0 - p1 = e1
            if net.Tpred > 1
                for tt = 2:net.Tpred
                    tempbat(1,tt) = xi.bus(1+Nabus*3+Ngen*2+2,tt-1) - xi.bus(1+Nabus*3+Ngen*2+1,tt) - xi.bus(1+Nabus*3+Ngen*2+2,tt);
                end
            end
            df.bus(1+Nabus*3+Ngen*2+1,:) = - dual{idx}.lambda(2+Nabus+1,:) - RHO * tempbat + dual{idx}.lambda(1,:) + RHO * temp1;
            df.bus(1+Nabus*3+Ngen*2+2,:) = - dual{idx}.lambda(2+Nabus+1,:) - RHO * tempbat + [dual{idx}.lambda(2+Nabus+1,2:net.Tpred) 0] + RHO * [tempbat(1,2:net.Tpred) 0];
        end
        
        % ancestor bus gradient 
        if Nabus > 0
            idxA = net.abus(idx,2);
            df.abus = dual{idx}.mu.abus + RHO * (xi.abus - z{idxA}(1,:)) ...
                - dual{idx}.lambda(3,:) - RHO * temp3; 
        end
        
        % children bus gradient 
        if Ncbus > 0
            for kC = 1:Ncbus
                idxC = net.cbus(idx,kC+1);
                df.cbus(1+3*(kC-1),:) = dual{idx}.mu.cbus(1+3*(kC-1),:) + RHO*(xi.cbus(1+3*(kC-1),:)-z{idxC}(2,:)) ...
                    - dual{idx}.lambda(1,:)*net.ri(idxC) - RHO*temp1*net.ri(idxC) ...
                    - dual{idx}.lambda(2,:)*net.xi(idxC) - RHO*temp2*net.xi(idxC);
                df.cbus(2+3*(kC-1),:) = dual{idx}.mu.cbus(2+3*(kC-1),:) + RHO*(xi.cbus(2+3*(kC-1),:)-z{idxC}(3,:)) ...
                    + dual{idx}.lambda(1,:) + RHO*temp1;
                df.cbus(3*kC,:) = dual{idx}.mu.cbus(3*kC,:) + RHO*(xi.cbus(3*kC,:)-z{idxC}(4,:)) ...
                    + dual{idx}.lambda(2,:) + RHO*temp2;
            end
        end
    end
end