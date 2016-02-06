function [Ai,bi,Di,Ei,Fi,LBi,UBi] = pre_opf_matrix(Ngen,Nabus,Ncbus,Tpred,Cst)

%:: xi = [v_i, L_i, P_i, Q_i, (p_i^gen,q_i^gen), (v_Ai), (L_Ci,P_Ci,Q_Ci)]

BD_INF = 1e3;

if Tpred == 1
    Nele = 4+Ngen*2+Nabus+Ncbus*3;
    Ai = zeros(3,Nele);
    bi = zeros(3,1);
    Ei = zeros(Nele,Nele);
    % objective 
    if size(Cst.cost,2) == 2
        Fi = zeros(Nele,1);
        if Ngen > 0 
            Fi(5,1) = Cst.cost(1,1);
        end
    else
        error('cost function order is larger than 1'); 
    end
    % constraint #1
    if Nabus > 0
        Ai(1,1) = 1;
        Ai(1,2) = Cst.ri^2 + Cst.xi^2;
        Ai(1,3) = - 2 * Cst.ri;
        Ai(1,4) = - 2 * Cst.xi;
        Ai(1,4+Ngen*2+1) = -1;
    end
    % constraint #2
    Ai(2,3) = - 1; 
    if Ngen > 0 
        Ai(2,5) = 1;
    end
    if Ncbus > 0
        for icbus = 1:Ncbus
            Ai(2,4+Ngen*2+Nabus+(icbus-1)*3+2) = 1;
            Ai(2,4+Ngen*2+Nabus+(icbus-1)*3+1) = - Cst.rc(icbus);
        end
    end
    % constraint #3
    Ai(3,4) = - 1;   
    if Ngen > 0 
        Ai(3,6) = 1;
    end
    if Ncbus > 0
        for icbus = 1:Ncbus
            Ai(3,4+Ngen*2+Nabus+(icbus-1)*3+3) = 1;
            Ai(3,4+Ngen*2+Nabus+(icbus-1)*3+1) = - Cst.xc(icbus);
        end
    end            
    % constraint #4
    if Nabus > 0
        Ei(1,2) = 1/2;
        Ei(2,1) = 1/2;
        Ei(3,3) = -1;
        Ei(4,4) = -1;
    end
    % consensus
    Di = eye(Nele);
    if Ngen > 0
        Di(5:6,:) = [];
    end
    % bounds
    LBi = - ones(Nele,1)*BD_INF;
    UBi = ones(Nele,1)*BD_INF;
    LBi(1,1) = Cst.bounds(1,1)^2;
    UBi(1,1) = Cst.bounds(1,2)^2;
    if Ngen > 0 
        LBi(5,1) = Cst.bounds(2,1);
        UBi(5,1) = Cst.bounds(2,2);
        LBi(6,1) = Cst.bounds(3,1);
        UBi(6,1) = Cst.bounds(3,2);
    end
else 
    error('Tpred > 1')
end

return

