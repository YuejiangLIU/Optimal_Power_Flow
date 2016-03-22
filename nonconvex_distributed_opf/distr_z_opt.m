function zi_new = distr_z_opt(idx,x,zi,dual,net,RHO)

% z_i update

zi_new = zeros(4,net.Tpred);
Nabus = net.abus(idx,1);
Ncbus = net.cbus(idx,1);

REGZ = 1e-10;        % REGZ * |zi_new - zi|^2

% minimize W*v^2 + U*v
W = RHO / 2 * (1+Ncbus) + REGZ;
U = - dual{idx}.mu.bus(1,:) - RHO * x{idx}.bus(1,:) - 2*REGZ*zi(1,:);
if Ncbus > 0 
    for kC = 1:Ncbus
        idxC = net.cbus(idx,kC+1);
        U = U - dual{idxC}.mu.abus - RHO * x{idxC}.abus;
    end
end
zi_new(1,:) = - U / (2*W);
clear W U

if Nabus > 0
    idxA = net.abus(idx,2);
    kC = find(net.cbus(idxA,:) == idx,1,'last') - 1;
    % minimize W*[L,P,Q]^2 + U*[L,P,Q]
    W = RHO + REGZ;
    U = - dual{idx}.mu.bus(2:4,:) - RHO * x{idx}.bus(2:4,:) ...
        - dual{idxA}.mu.cbus(1+(kC-1)*3:kC*3,:) - RHO * x{idxA}.cbus(1+(kC-1)*3:kC*3,:) ...
        - 2*REGZ*zi(2:4,:) ;
    zi_new(2:4,:) = - U / (2*W);
else 
    zi_new(2:4,:) = zeros(3,net.Tpred);
end

end