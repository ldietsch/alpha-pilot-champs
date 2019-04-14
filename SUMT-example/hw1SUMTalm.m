function A = hw1SUMTalm(x,r_p,lambda,c)
% This function is the pseudo-objective function using the ALM.
% In this function, r_p and lambda are parameters, and x are the 
% variables.
% Prof. Crossley 20 Sep 2016
% Updated by Luis Dietsche 2 Oct 2016

% compute values of the objective function and constraints at the current
% value of x
f = hw1SUMTfun(x);
g = hw1SUMTcon(x);
l = length(g);

L = 0;
psi = zeros(l,1);
% Fletcher's substitution
for i = 1:l
%     psi(i) = max(c(i)*g(i), -lambda(i) / (2 * r_p));
    psi(i) = max(g(i), -lambda(i) / (2 * r_p));
    L = L + lambda(i)*psi(i)+r_p*psi(i)^2;
end

 A = f+L;
end