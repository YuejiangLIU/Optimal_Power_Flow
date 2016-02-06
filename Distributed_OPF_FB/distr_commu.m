function zi_hat = distr_x_commu(abus,cbus,z,i,rhoc)

z_Ai = [];
z_Ci = [];

if abus(1) > 0
    idxAbus = abus(2);
    z_Ai= z{idxAbus}(1);
end
if cbus(1) > 0
    idxCbus = cbus(2:end);
    z_Ci = z{idxCbus}(2:4);
    if cbus(1) > 1
        error('check z_Ci');
    end
end

zi_hat = [z{i};z_Ai;z_Ci];

end