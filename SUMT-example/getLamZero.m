function lam0_eq = getLamZero(x0, U, f, h, info)
% getLamZero sets initial values for Lagrange multipliers related to the
% SUMT method of Augmented Lagrange Multipliers
% f and h are function handles for the Lagrangian and equality constraints,
% respectively
% delf is a (N*dimM) x 1 vector, where N is the number of segments in the
% trajectory and dimM is the dimension of the control input in state space
% form
delf = getPartialDerivs(x0, U, f, info);
% delh is a k x (N*dimM) vector, where k is the number of equality
% constraints
delh = getPartialDerivs(x0, U, h, info);

k = info.k_eq; % number of equality constraints
lam0_eq = zeros(k,1);
for i = 1:k
    s = dot(delf, delh(i,:));
    if s <= 0
        lam0_eq(i) = 1;
    else
        lam0_eq(i) = -1;
    end
end

end

