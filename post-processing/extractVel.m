function vel = extractVel(DX, N)

vel = zeros(3*N,1);
j = 1;
k = 1;
for i = 1:N
    vel(j:3*i) = DX(k:k+2);
    j = j+3;
    k = k + 12;
end

vel = reshape(vel, 3, N);

end
