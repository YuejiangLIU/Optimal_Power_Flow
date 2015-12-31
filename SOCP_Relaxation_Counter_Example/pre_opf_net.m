function net = pre_opf_net(mpc)

% parameter
net.Nbus = size(mpc.bus,1);
net.Nline = size(mpc.branch,1);

% demand
net.p_demand = mpc.bus(:,3)/mpc.baseMVA;
net.q_demand = mpc.bus(:,4)/mpc.baseMVA;

% generator 
net.p_gen_max = zeros(net.Nbus,1);
net.p_gen_min = zeros(net.Nbus,1);
net.q_gen_max = zeros(net.Nbus,1);
net.q_gen_min = zeros(net.Nbus,1);

net.ncost = mpc.gencost(1,4);
net.c_gen = zeros(net.Nbus,net.ncost);

for ii = 1:size(mpc.gen,1)
    idxGen = mpc.gen(ii,1);
    net.p_gen_max(idxGen) = mpc.gen(ii,9)/mpc.baseMVA;
    net.p_gen_min(idxGen) = mpc.gen(ii,10)/mpc.baseMVA;
    net.q_gen_max(idxGen) = mpc.gen(ii,4)/mpc.baseMVA;
    net.q_gen_min(idxGen) = mpc.gen(ii,5)/mpc.baseMVA;
    net.c_gen(idxGen,:) = mpc.gencost(ii,end-net.ncost+1:end);
end

% shunt 
net.G_shunt = mpc.bus(:,5)/mpc.baseMVA;
net.B_shunt = mpc.bus(:,6)/mpc.baseMVA;

% impedance
net.r = mpc.branch(:,3);
net.x = mpc.branch(:,4);

% nodal admittance
net.G = zeros(net.Nbus,net.Nbus);
net.B = zeros(net.Nbus,net.Nbus);

for ii = 1:net.Nline
    idxFrom = mpc.branch(ii,1);
    idxTo = mpc.branch(ii,2);
    z = complex(mpc.branch(ii,3),mpc.branch(ii,4));
    y = 1/z;
    net.G(idxFrom,idxTo) = - real(y);
    net.B(idxFrom,idxTo) = - imag(y);
    net.G(idxTo,idxFrom) = - real(y);
    net.B(idxTo,idxFrom) = - imag(y);     
end

for ii = 1:net.Nbus 
    net.G(ii,ii) = net.G_shunt(ii) - sum(net.G(ii,:));
    net.B(ii,ii) = net.B_shunt(ii) - sum(net.B(ii,:));
end

% voltage
net.v_max = mpc.bus(:,12);
net.v_min = mpc.bus(:,13);

idxRef = mpc.bus(:,2) == 3 ;                % referance bus
net.v_max(idxRef) = mpc.bus(idxRef,8);
net.v_min(idxRef) = mpc.bus(idxRef,8);

end