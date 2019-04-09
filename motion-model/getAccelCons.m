function accel_cons = getAccelCons(acc, a_max, N)

accel_cons = zeros(N,1);
for i = 1:N
   
    accel_cons(i) = (norm(acc(:,i)) - a_max)/a_max;
    
end


end