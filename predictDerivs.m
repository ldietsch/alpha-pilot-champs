function DX = predictDerivs(x, u, N, info)

DX = zeros(12*N,1);
j = 1;
k = 1;
for i = 1:N
   DX(j:i*12) = getStateDerivs(x(j:i*12), u(k:i*4), info);
   j = j + 12;
   k = k + 4;
end

end