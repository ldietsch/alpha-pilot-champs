function pos = extractPosVec(results, steps)

pos = zeros(3*(steps+1),1);
j = 1;
k = 1;
for i = 1:steps+1
    pos(j:3*i) = results(k:k+2);
    j = j+3;
    k = k + 12;
end

end