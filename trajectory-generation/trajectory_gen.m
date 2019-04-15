x0 = [0;-1;0];    % initial design
c = findc(x0);
iter = 1;         % initial value of minimization counter 
rp = 1.0;     % initial value of penalty multiplier 
g = hw1SUMTcon(x0);
l = length(g);
lambda = zeros(l,1);     % initial Lagrange multipliers
gam = 3;
% compute function value at x0, initialize convergence criteria
fnew = hw1SUMTfun(x0);
fold = 2 * fnew;   % ensure that first loop does not trigger convergence
% set optimization options - use default BFGS with numerical gradients
% provide display each iteration
options = optimoptions(@fminunc,'Algorithm', 'quasi-newton', ...
 'Display', 'iter');

% begin sequential minimizations  - note tolerances chosen here
% absolute tolerance for change in objective function (1e-3, here) and 
% absolute tolerance for constraints (1e-5, here)
maxiter = 100;
t3  = table();
options.StepTolerance = 1e-9 ;
while abs(fnew-fold) >= 1e-6 || (max(g) >= 1e-6) && iter<maxiter
    fold = fnew;  % store last objective function value
    % call fminunc - use "ALM" pseudo-objective function, note that r_p and
    % lambda are passed as parameters, no semi-colon to display results
    [xstar,phistar,exitflag] = fminunc(@(x)hw1SUMTalm(x,rp,lambda,c),x0,options);
    % compute objective and constraints at current xstar
    f = hw1SUMTfun(xstar);
    fnew = f;
    g = hw1SUMTcon(xstar);
    % update lagrange multipliers
    for i = 1:l
        lambda(i) = lambda(i) + 2 * rp * max( g(i), -lambda(i)/(2*rp));
    end
    t3 = [t3;table(rp,x0',xstar',f,g',iter,exitflag)];
    iter = iter + 1;     % increment minimization counter
    rp = rp * gam; % increase penalty multiplier (gamma = 5, here)
    x0 = xstar;    % use current xstar as next x0
    
end