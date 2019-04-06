function bAngAcc = getBodyAngAcc(bAngV, tau, J)
tau_phi = tau(1); tau_theta = tau(2); tau_psi = tau(3);
Jx = J(1); Jy = J(2); Jz = J(3);

bAngAcc = [(Jy-Jz)*q*r/Jx; (Jz - Jx)*p*r/Jy; (Jx - Jy)*p*q/Jz] + ...
    [tau_phi/Jx; tau_theta/Jy; tau_psi/Jz];

end