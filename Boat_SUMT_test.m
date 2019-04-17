clear, clc, close all
format long    % use format long to see differences near convergence

N = 20;

U0 = ones(2,N);% initial design
Ts = [0.05; 0];
U0 = [U0, Ts];

iter = 1;         % initial value of minimization counter 
rp = 1.0;     % initial value of penalty multiplier 

g = inequalCon(U0, N);
X = intBoat(U0, N);
h = equalCon(X, N);

l = length(g);
lambda_g = zeros(4,N+1);     % initial Lagrange multipliers
lambda_h = [-1 ; 1];
gam = 3;
% compute function value at x0, initialize convergence criteria
fnew = costFun(U0, N);
fold = 2 * fnew;   % ensure that first loop does not trigger convergence
% set optimization options - use default BFGS with numerical gradients
% provide display each iteration
options = optimoptions(@fminunc,'Algorithm', 'quasi-newton', ...
 'Display', 'iter');

% begin sequential minimizations  - note tolerances chosen here
% absolute tolerance for change in objective function (1e-3, here) and 
% absolute tolerance for constraints (1e-5, here)
maxiter = 100;

while abs(fnew-fold) >= 1e-9 && iter<maxiter
    fold = fnew;  % store last objective function value
    % call fminunc - use "ALM" pseudo-objective function, note that r_p and
    % lambda are passed as parameters, no semi-colon to display results
    [Ustar,phistar,exitflag] = fminunc(@(U)psuedoFun(U, rp, lambda_g, lambda_h, N),U0,options);
    % compute objective and constraints at current xstar
    f = costFun(Ustar, N);
    fnew = f;
    g = inequalCon(Ustar, N);
    X = intBoat(Ustar, N);
    h = equalCon(X, N);


    % update lagrange multipliers
    for j = 1:4
        for i = 1:N+1
            lambda_g(j,i) = lambda_g(j,i) + 2 * rp * max( g(j,i), -lambda_g(j,i)/(2*rp));
        end
        
    end
    
    for j = 1:2
        lambda_h(j) = lambda_h(j) + 2*rp*h(j);
    end
    iter = iter + 1;     % increment minimization counter
    rp = rp * gam; % increase penalty multiplier (gamma = 5, here)
    U0 = Ustar;    % use current xstar as next x0
    
end
%X = intBoat(Ustar, N);


x = X(1,:); y = X(2,:);
plot(x,y)

t = 0:Ts:(N-1)*Ts;
figure
%plot(t,(Ustar(2,1:N))*180/pi)
plot(t,y)

Ustar(1,N+1)
