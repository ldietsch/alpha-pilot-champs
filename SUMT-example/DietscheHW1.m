 clear, clc
%Luis Dietsche
%Multi-Disciplinary Design Optimization
%10/15/2018
%% Part I: Engineering Problems in N Variables
% This part of the homework utilizes various optimization techniques to
% obtain the minimization of a potential energy function of a three-bar
% truss problem using two dofs in horizontal and vertical
% displacements
x0 = [0;0]; %initial guess the same for all options
%% Question 3
%Solution with finite difference gradients
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient', false, 'Display', 'iter');
[x1, fval1, exitflag1, output1] = fminunc(@PiU, x0, options);
[f1, del1] = PiU(x1);
%solution with analytic gradients
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient', true, 'Display', 'iter');
[x2, fval2, exitflag2, output2] = fminunc(@PiU, x0, options);
[f2, del2] = PiU(x2);
%% Question 4
%Solution with analytic gradients and DFP update
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
'SpecifyObjectiveGradient', true, 'Display', ...
'iter', 'HessUpdate','dfp');
[x3, fval3, exitflag3, output3] = fminunc(@PiU, x0, options);
[f3, del3] = PiU(x3);
%Solution with analytic gradients and steepest descent update
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
'SpecifyObjectiveGradient', true, 'Display', ...
'iter', 'HessUpdate','steepdesc');
[x4, fval4, exitflag4, output4] = fminunc(@PiU, x0, options);
[f4, del4] = PiU(x4);
%% Question 5
%Modified Newton's Method with analytic gradient and Hessian
options = optimoptions(@fminunc, 'Algorithm',...
'trust-region', 'SpecifyObjectiveGradient', true, 'Display', 'iter',...
'HessianFcn', 'objective');
[x5, fval5, exitflag5, output5] = fminunc(@PiU, x0, options);
[f5, del5] = PiU(x5);
pause

%% Part II: Engineering Application of SUMT Approach
% DESCRIPTIVE TEXT
clear, clc
%write feasible x0 (necessary for interior penalty), keeping in mind to
%check lift and lift coefficient constraints
x0 = [2;12;8*pi/180]; 
%% 2a Exterior Penalty Method
%find scaling factors cj =  del(f(x0))/del(gj(x0)), do cheap approximation
c = findc(x0);

fold = hw1SUMTphi(x0,1,c);
%get x1
%Solution with finite difference gradients
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient', false);
[x1, fnew, exitflag, output] = fminunc(@(x)hw1SUMTphi(x,1,c),x0, options );
g = hw1SUMTcon(x1);
f = hw1SUMTfun(x1);
rp = 1;
xold = x1;
eps = 10^-4;
rp = 1; gam = 3;
maxiter = 100;
iter = 0;
options.StepTolerance = 1e-9;
%Exterior penalty method with scaling of constraints
t0 = table(rp,x0',x1',f,g',iter,exitflag);
while abs(fold-fnew)>=eps && iter<=maxiter
    fold = fnew;
    rp = gam*rp;
    [xnew, fnew, exitflag, output] = ...
        fminunc(@(x)hw1SUMTphi(x,rp,c),xold,options );
    g = hw1SUMTcon(xold);
    f = hw1SUMTfun(xold);
    iter = iter+1;
    t0 = [t0;table(rp,xold',xnew',f,g',iter,exitflag)];
    xold = xnew;
end
format long
t0
[D L CL] = checkxstar(xnew)
pause
%% 2b Interior Penalty Method
clearvars -except t0 c, clc
x0 = [2;12;8*pi/180]; 
c = findc(x0);
fold = hw1SUMTphiInt(x0,1,c);
%get x1
%Solution with finite difference gradients
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient',false,'Display','iter');
[x1, fnew, exitflag, output] = fminunc(@(x)hw1SUMTphiInt(x,1,c),x0, options );
g = hw1SUMTcon(x1);
f = hw1SUMTfun(x1);
xold = x1;
eps = 10^-4;
rp = 1000; gam = 10;
maxiter = 100;
iter = 1;
options.StepTolerance = 1e-9 ;
%Interior penalty method with scaling of constraints
t1 = table(rp,x0',x1',f,g',iter,exitflag);
while abs(fold-fnew)>=eps && iter<=maxiter
    fold = fnew;
    rp = (1/gam)*rp;
    [xnew, fnew, exitflag, output] = ...
        fminunc(@(x)hw1SUMTphiInt(x,rp,c),xold,options );
    g = hw1SUMTcon(xold);
    f = hw1SUMTfun(xold);
    iter = iter+1;
    t1 = [t1;table(rp,xold',xnew',f,g',iter,exitflag)];
    xold = xnew;

end   
    t1
    [D L CL] = checkxstar(xnew)
    pause
%% 2c Linear Extended Method    
clearvars -except t0 t1, clc
x0 = [0;-1;0]; 
eps0 = -0.2;
rp = findrp(x0,eps0);
C = (-eps0)/rp^1/2;
fold = hw1SUMTfun(x0);

c = findc(x0);

%get x1
%Solution with finite difference gradients
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient',false,'Display','iter');
[x1, fnew, exitflag, output] = fminunc(@(x)hw1SUMTphiLinExt(x,rp,eps0),x0, options );
g = hw1SUMTcon(x1);
f = hw1SUMTfun(x1);
xold = x1;
eps = 10^-4;
gam = 4;
maxiter = 100;
iter = 1;
options.StepTolerance = 1e-9 ;
%Interior penalty method with scaling of constraints
t2 = table(rp,x0',x1',f,g',iter,exitflag);
while abs(fold-fnew)>=eps && iter<=maxiter
    fold = fnew;
    rp = (1/gam)*rp;
    eps0 = -C*rp^1/2;
    [xnew, fnew, exitflag, output] = ...
        fminunc(@(x)hw1SUMTphiLinExt(x,rp,eps0),xold,options );
    g = hw1SUMTcon(xold);
    f = hw1SUMTfun(xold);
    iter = iter+1;
    t2 = [t2;table(rp,xold',xnew',f,g',iter,exitflag)];
    xold = xnew;

end   
    t2
[D L CL] = checkxstar(xnew)
pause
%% Problem 2d Augmented Lagrange Multipliers
clearvars -except t0 t1 t2, clc
format long    % use format long to see differences near convergence

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
t3
[D L CL] = checkxstar(xstar)
%% Functions    
function rp = findrp(x,eps0)
f = hw1SUMTfun(x);
g = hw1SUMTcon(x);
    %linear extended method
ncon = length(g);   % number of constraints
P = 0;              % intialize P value to zero
for j = 1:ncon
    if g(j) <= eps0
        g(j) = -1/(g(j));
    else
        g(j) = -(2*eps0-g(j))/eps0^2;
    end
      P = P + g(j);
end
rp = abs(f)/P;

end
    
function c = findc(x0)
    
delta = 10^-4;
g0 = hw1SUMTcon(x0);
g0delt = hw1SUMTcon(x0+delta);
f0 = hw1SUMTfun(x0);
f0delt = hw1SUMTfun(x0+delta);
n = length(g0);
c = zeros(n,1);
for i = 1:n
    c(i) = abs(f0delt-f0)/abs(g0delt(i)-g0(i));
end   
    
end
    
    
    


