function vel_cons = getVelCons(vel, v_max, N)

vel_cons = zeros(3,2*N);
j = N+1;
for i = 1:N
   
    vel_cons(:,i) = (vel(:,i) - v_max)./v_max;
    vel_cons(:,j) = (-vel(:,i) - v_max)./v_max;
    j = j+1;
end


end