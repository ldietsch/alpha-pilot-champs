function fx = getStateDerivs(x,u,info)
%getStates(x,u) takes the current states and inputs and computes the first 
%order derivatives of each state. The order of the output is [pxd pyd pzd 
%udb vdb wdb phid thetad psid pdb qdb rdb]. The order of the current states
%are [px py pz ub vb wb phi theta psi pb qb rb]. Motor inputs are in the
%order of front, right, back, left.
m = info.m;
g = info.g;
J = info.J;
l = info.l;
k_tau = info.k_tau;
F = sum(u);
%p_inert = x(1:3); % not used here
v_body = x(4:6);
eAngles = x(7:9);
bAngV = x(10:12);
v_inert = getVelInert(eAngles, v_body);
acc_body = getAccBody(bAngV, v_body, eAngles, m, F, g);
vEangles = getVelEangles(eAngles, bAngV);
tau = getBodyTorque(u, l, k_tau); 
bAngAcc = getBodyAngAcc(bAngV, tau, J);

fx = [v_inert; acc_body; vEangles; bAngAcc];


end

