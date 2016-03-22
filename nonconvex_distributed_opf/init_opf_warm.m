function [x_init,z_init,dual_init] = init_opf_warm(x,z,dual,net) 

% :: initialize the new algorithm from the old solution by shifting one hour

for idx = 1:net.Nbus
    Nabus = net.abus(idx,1);
    Ncbus = net.cbus(idx,1);
    Ngen = net.gen(idx,1);
    Nbat = net.bat(idx,1);
    
    x_init{idx}.bus = circshift(x{idx}.bus,-1,2);
    z_init{idx} = circshift(z{idx},-1,2);
    dual_init{idx}.lambda = circshift(dual{idx}.lambda,-1,2);
    dual_init{idx}.mu.bus = circshift(dual{idx}.mu.bus,-1,2);
    
    if Nabus > 0 
        x_init{idx}.abus = circshift(x{idx}.abus,-1,2);
        dual_init{idx}.mu.abus = circshift(dual{idx}.mu.abus,-1,2);
        dual_init{idx}.gamma = circshift(dual{idx}.gamma,-1,2);
    end
    
    if Ncbus > 0 
        x_init{idx}.cbus = circshift(x{idx}.cbus,-1,2);
        dual_init{idx}.mu.cbus = circshift(dual{idx}.mu.cbus,-1,2);
    end        
    
end

end