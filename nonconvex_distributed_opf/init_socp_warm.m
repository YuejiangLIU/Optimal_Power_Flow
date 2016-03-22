function [x_init,z_init,dual_init] = init_socp_warm(sol_ref,net)

% :: Intialization from convex relaxed opf solution 

for idx = 1:net.Nbus

    Nabus = net.abus(idx,1);
    Ncbus = net.cbus(idx,1);
    Ngen = net.gen(idx,1);
    Nbat = net.bat(idx,1);
    
    x_init{idx}.bus(1,:) = sol_ref.v(idx,:);
    z_init{idx}(1,:) = sol_ref.v(idx,:);
    
    if Nabus > 0
        x_init{idx}.bus(2,:) = sol_ref.L(idx,:);
        x_init{idx}.bus(3,:) = sol_ref.P(idx,:);
        x_init{idx}.bus(4,:) = sol_ref.Q(idx,:);
    end
    z_init{idx}(2,:) = sol_ref.L(idx,:);
    z_init{idx}(3,:) = sol_ref.P(idx,:);
    z_init{idx}(4,:) = sol_ref.Q(idx,:);

    if Ngen > 0
        x_init{idx}.bus(1+Nabus*3+1,:) = sol_ref.p_gen(idx,:);
        x_init{idx}.bus(1+Nabus*3+2,:) = sol_ref.q_gen(idx,:);
    end
    
    if Nbat > 0 
        x_init{idx}.bus(1+Nabus*3+Ngen*2+1,:) = sol_ref.p_bat(idx,:);
        x_init{idx}.bus(1+Nabus*3+Ngen*2+2,:) = sol_ref.e_bat(idx,:);
    end
    
    if Nabus > 0
        idxA = net.abus(idx,2);
        x_init{idx}.abus = sol_ref.v(idxA,:);
    end

    if Ncbus > 0
        for kC = 1:Ncbus
            idxC = net.cbus(idx,kC+1);
            x_init{idx}.cbus(1+(kC-1)*3,:) = sol_ref.L(idxC,:);
            x_init{idx}.cbus(2+(kC-1)*3,:) = sol_ref.P(idxC,:);
            x_init{idx}.cbus(kC*3,:) = sol_ref.Q(idxC,:);
        end
    end
    
    % Dual Initialization
    dual_init{idx}.lambda = zeros(2+Nabus+Nbat,net.Tpred);
    dual_init{idx}.gamma = zeros(Nabus,net.Tpred);
    dual_init{idx}.mu.bus = zeros(1+Nabus*3,net.Tpred);
    if Nabus > 0
        dual_init{idx}.mu.abus = zeros(Nabus,net.Tpred);
    end
    if Ncbus > 0
        dual_init{idx}.mu.cbus = zeros(Ncbus*3,net.Tpred);
    end
end
    
end