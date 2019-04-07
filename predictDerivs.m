function DX = predictDerivs(x, u, N, info)

DX = zeros(12*N,1);
j = 1;
for i = 1:N
   DX(j:i*12) = getStateDerivs(x(j:i*12), u(j:i*12), info);
   j = j + 12;
end

end