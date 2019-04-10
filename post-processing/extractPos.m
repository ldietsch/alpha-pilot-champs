function [x_i, y_i, z_i] = extractPos(results, steps)

pos = zeros(3*steps,1);
j = 1;
k = 1;
for i = 1:steps
    pos(j:3*i) = results(k:k+2);
    j = j+3;
    k = k + 12;
end

pos = reshape(pos, 3, steps);
x = pos(1,:)
end