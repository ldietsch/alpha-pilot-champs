function phi = hw1SUMTphi(x,r_p,c)
% This function is the pseudo-objective function using the exterior penalty.
% In this function, r_p is a "parameter", x are the variables. 
% Prof. Crossley 20 Sep 2016

% compute values of the objective function and constraints at the current
% value of x
f = hw1SUMTfun(x);
g = hw1SUMTcon(x);

% exterior penalty function
ncon = length(g);   % number of constraints
P = 0;              % intialize P value to zero
for j = 1:ncon
    P = P + c(j)*max(0,g(j))^2; 
%     P = P + max(0,g(j))^2; 
end
phi = f + r_p * P;

end