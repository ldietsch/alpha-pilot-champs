function phi = hw1SUMTphiLinExt(x,r_p,eps)
% This function is the pseudo-objective function using the linear extended
% interior penalty function.
% In this function, r_p is a "parameter", x are the variables. 

% compute values of the objective function and constraints at the current
% value of x
f = hw1SUMTfun(x);
g = hw1SUMTcon(x);

% linear extended interior penalty function
ncon = length(g);   % number of constraints
P = 0;              % intialize P value to zero
for j = 1:ncon
    if g(j) <= eps
        g(j) = -1/(g(j));
    else

        g(j) = -(2*eps-g(j))/eps^2;
    end
      P = P + g(j);
end
phi = f + r_p * P;

end