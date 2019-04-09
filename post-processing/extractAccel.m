function acc = extractAccel(DX, N)

acc = zeros(3*N,1);
j = 1;
k = 4;
for i = 1:N
    acc(j:3*i) = DX(k:k+2);
    j = j+3;
    k = k + 12;
end

acc = reshape(acc, 3, N);

end

