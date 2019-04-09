function Qbar = getQbar(Q0, n, N)
Qbar = zeros(n*N, n*N);
j = 1;
k = n;
for i = 1:N
    Qbar(j:k, j:k) = blkdiag(Q0);
    j = j + n;
    k = k + n;
end

end