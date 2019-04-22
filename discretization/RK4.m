function fx_k_1 = RK4(x, u, Ts, info)


k = zeros(12,4);
fx_k1 = getStateDerivs(x,u,info);
k(:,1) = Ts*fx_k1;
fx_k2 = getStateDerivs(x+k(:,1)/2,u,info);
k(:,2) = Ts*fx_k2;
fx_k3 = getStateDerivs(x+k(:,2)/2,u,info);
k(:,3) = Ts*fx_k3;
fx_k4 = getStateDerivs(x+k(:,3), u, info);
k(:,4) = Ts*fx_k4;
fx_k_1 =  x + 1/6*(k(:,1) + 2*k(:,2) + 2*k(:,3) + k(:,4));
    
end