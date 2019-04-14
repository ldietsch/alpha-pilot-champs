function phi = hw1SUMTphiInt(x,r_p,c)
% This function is the pseudo-objective function using the exterior penalty.
% In this function, r_p is a "parameter", x are the variables. 

% compute values of the objective function and constraints at the current
% value of x
f = hw1SUMTfun(x);
g = hw1SUMTcon(x);

% exterior penalty function
ncon = length(g);   % number of constraints
P = 0;              % intialize P value to zero
for j = 1:ncon
%     P = P + (-log(-g(j)));
      P = P + (-log(-c(j)*g(j)));
end
phi = f + r_p * P;

end