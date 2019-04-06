function tau = getBodyTorque(u, l, k_tau)
Ff = u(1); Fr = u(2); Fb = u(3); Fl = u(4);
tau_phi = l*(Fl - Fr);
tau_theta = l*(Ff - Fb);
tau_psi = k_tau*(Fr+Fl-Ff-Fb);

tau = [tau_phi;tau_theta;tau_psi];

end