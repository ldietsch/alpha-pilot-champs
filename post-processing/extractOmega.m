function omega = extractOmega(results, steps)

omega = zeros(3*steps,1);
k = 10;
j =1;
for i = 1:steps
    omega(j:3*i) = results(k:k+2);
    j=j+3;
    k = k + 12;
end

end