function [di_new,res_sum_square] = distr_dual_update(idx,xi,z,di,net,par,RHO)

% :: dual varaibles update

Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);
Ngen = net.gen(idx,1);
Nbat = net.bat(idx,1);

temp1 = - par.p_demand(idx,:);
temp2 = - par.q_demand(idx,:);
if Nabus > 0
    temp1 = temp1 - xi.bus(3,:);
    temp2 = temp2 - xi.bus(4,:);
    temp3 = xi.bus(1,:) + (net.ri(idx)^2+net.xi(idx)^2)*xi.bus(2,:) - 2*(net.ri(idx)*xi.bus(3,:)+net.xi(idx)*xi.bus(4,:)) - xi.abus;
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

di_new.lambda(1:2,:) = di.lambda(1:2,:) + RHO * [temp1; temp2];
if Nabus > 0 
    di_new.lambda(3,:) = di.lambda(3,:) + RHO * temp3;
    di_new.gamma = di.gamma + RHO * tempnc;
end
if Nbat > 0
    tempbat(1,1) = net.e_bat_now(idx) - xi.bus(1+Nabus*3+Ngen*2+1,1) - xi.bus(1+Nabus*3+Ngen*2+2,1);    % e0 - p1 = e1
    for tt = 2:net.Tpred
        tempbat(1,tt) = xi.bus(1+Nabus*3+Ngen*2+2,tt-1) - xi.bus(1+Nabus*3+Ngen*2+1,tt) - xi.bus(1+Nabus*3+Ngen*2+2,tt);
    end
    di_new.lambda(2+Nabus+1,:) = di.lambda(2+Nabus+1,:) + RHO * tempbat;
end
if Nabus > 0
    diffbus = xi.bus(1:4,:) - z{idx}(1:4,:);
else
    diffbus = xi.bus(1,:) - z{idx}(1,:);
end
di_new.mu.bus = di.mu.bus + RHO * diffbus;
if Nabus > 0
    idxA = net.abus(idx,2);
    diffabus = xi.abus -z{idxA}(1,:);
    di_new.mu.abus = di.mu.abus + RHO * diffabus;
end
if Ncbus > 0
    for kC = 1:Ncbus
        idxC = net.cbus(idx,kC+1);
        diffcbus(1+(kC-1)*3:kC*3,:) = xi.cbus(1+(kC-1)*3:kC*3,:) - z{idxC}(2:4,:);
    end
    di_new.mu.cbus = di.mu.cbus + RHO * diffcbus;
end

% residual 
if Nabus > 0 
    res_sum_square(1,1) = temp1*temp1' + temp2*temp2' + temp3*temp3' + tempnc*tempnc';
else
    res_sum_square(1,1) = temp1*temp1' + temp2*temp2';
end
if Nbat > 0 
    res_sum_square(1,1) = res_sum_square(1,1) + tempbat*tempbat';
end
res_sum_square(1,2) = sum(sum(diffbus.^2));
if Nabus > 0 
    res_sum_square(1,2) = res_sum_square(1,2) + diffabus*diffabus';
end
if Ncbus > 0 
    res_sum_square(1,2) = res_sum_square(1,2) + sum(sum(diffcbus.^2));
end

end