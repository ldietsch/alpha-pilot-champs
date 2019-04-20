function angles = extractAngles(X, N)

angles = zeros(3*N,1);
j = 1;
k = 7;
for i = 1:N
    angles(j:3*i) = X(k:k+2);
    j = j+3;
    k = k + 12;
end

angles = reshape(angles, 3, N);

end