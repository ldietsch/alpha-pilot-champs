function vel_cons = getVelCons(vel, v_max, N)

vel_cons = zeros(3,N);
for i = 1:N
   
    vel_cons(:,i) = (vel(:,i) - v_max)./v_max;
    
end


end