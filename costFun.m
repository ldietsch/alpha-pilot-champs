function f = costFun(U, N)

rho = .001;
cost_u = 0;
Ts = U(1,N+1);
for i = 1:N
   cost_u = cost_u + U(:,i)'*U(:,i)*rho;  
end
f = cost_u + N*Ts^2;