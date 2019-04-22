function psi_euler = extractHeading(results, steps)

psi_euler = zeros(steps,1);
k = 1;
for i = 1:steps
    psi_euler(i) = results(k+8);
    k = k + 12;
end

end