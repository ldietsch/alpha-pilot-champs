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