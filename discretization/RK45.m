function [fx_k_1, Ts_new] = RK45(x, u, info)
Ts = info.Ts;
tol = info.tol;
k = zeros(12,6);
fx_k1 = getStateDerivs(x,u,info);
k(:,1) = Ts*fx_k1;
fx_k2 = getStateDerivs(x+k(:,1)/4,u,info);
k(:,2) = Ts*fx_k2;
fx_k3 = getStateDerivs(x+k(:,1)*3/32 + 9*k(:,2)/32,u,info);
k(:,3) = Ts*fx_k3;
fx_k4 = getStateDerivs(x+1932*k(:,1)/2197 - 7200*k(:,2)/2197 + ...
    7296*k(:,3)/2197, u, info);
k(:,4) = Ts*fx_k4;
fx_k5 = getStateDerivs(x + 439*k(:,1)/216 - 8*k(:,2) + ...
    3680*k(:,3)/513 - 845*k(:,4)/4104, u, info);
k(:,5) = Ts*fx_k5;
fx_k6 = getStateDerivs(x - 8*k(:,1)/27 + 2*k(:,2) - 3544*k(:,3)/2565 + ...
    1859*k(:,4)/4104 - 11*k(:,5)/40, u, info);
k(:,6) = Ts*fx_k6;
xrk4 = x + 25*k(:,1)/216 + 1408*k(:,3)/2565 + 2197*k(:,4)/4101 - ...
    1/5*k(:,5);
fx_k_1 =  x + 16/135*k(:,1) + 6656/12825*k(:,3) + 28561/56430*k(:,4) - ...
    9/50*k(:,5) + 2/55*k(:,6);
s = (tol/(2*norm(fx_k_1 -  xrk4)))^(1/4);
Ts_new = s*Ts;

end