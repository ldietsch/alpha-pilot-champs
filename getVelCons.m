function vel_cons = getVelCons(vel, v_max, N)

vel_cons = zeros(N,1);
for i = 1:N
   
    vel_cons(i) = (norm(vel(:,i)) - v_max)/v_max;
    
end


end