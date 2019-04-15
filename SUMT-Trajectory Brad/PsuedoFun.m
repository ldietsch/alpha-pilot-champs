function phi = PsuedoFun(f, h)

f = MinFunc(x);
g = inequalCon(x);

% exterior penalty function
ncon = length(g);   % number of constraints
P = 0;              % intialize P value to zero
for j = 1:ncon
    P = P + c(j)*max(0,g(j))^2; 
%     P = P + max(0,g(j))^2; 
end
phi = f + r_p * P;