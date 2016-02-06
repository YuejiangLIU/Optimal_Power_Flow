function zi = distr_z_opt(x,abus,cbus,dual,i,Ngen,RHO)

zi = zeros(4,1);
Nabus = abus(i,1);
Ncbus = cbus(i,1);

if Ncbus < 1
    zi(1,1) = x{i}(1) + dual{i}.mu(1)/RHO;
else     
    % min W*v^2 + U*v
    W = RHO / 2 * (1+Ncbus);
    idxcbus = cbus(i,2:Ncbus+1); 
    elecbus = 5;
    U = - dual{i}.mu(1) - RHO * x{i}(1) ...
        - dual{idxcbus}.mu(elecbus) - RHO * x{idxcbus}(elecbus + Ngen(idxcbus)*2);
    zi(1,1) = - U/(2*W);
end 

clear W U
if Nabus < 1
    zi(2:4,1) = x{i}(2:4) + dual{i}.mu(2:4)/RHO;
else
    % min W*v^2 + U*v
    W = RHO;
    idxabus = abus(i,2);
    eleabus = [5:7]+abus(idxabus,1);
    U = - dual{i}.mu(2:4) - RHO * x{i}(2:4) ...
        - dual{idxabus}.mu(eleabus) - RHO * x{idxabus}(eleabus+Ngen(idxabus)*2);
    zi(2:4,1) = - U/(2*W);
end

end
