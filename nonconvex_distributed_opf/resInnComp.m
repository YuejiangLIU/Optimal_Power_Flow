function [x_res,z_res] = resInnComp(x_new,z_new,x,z,net)

% :: Inner loop residual x_res = |x_new - x|, z_res = |z_new - z|

x_res_sum_square = 0;
z_res_sum_square = 0;

for idx = 1:net.Nbus
    Nabus = net.abus(idx,1);
    Ncbus = net.cbus(idx,1);
    x_res_sum_square = x_res_sum_square + sum(sum( (x_new{idx}.bus - x{idx}.bus).*(x_new{idx}.bus - x{idx}.bus) )); 
    z_res_sum_square = z_res_sum_square + sum(sum( (z_new{idx} - z{idx}).*(z_new{idx} - z{idx}) ));
    if Nabus > 0  
        x_res_sum_square = x_res_sum_square + sum(sum( (x_new{idx}.abus - x{idx}.abus).*(x_new{idx}.abus - x{idx}.abus) ));
    end
    if Ncbus > 0 
        x_res_sum_square = x_res_sum_square + sum(sum( (x_new{idx}.cbus - x{idx}.cbus).*(x_new{idx}.cbus - x{idx}.cbus) ));
    end
end

x_res = sqrt(x_res_sum_square);
z_res = sqrt(z_res_sum_square);


end