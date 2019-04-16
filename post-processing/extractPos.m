function [x_i, y_i, z_i] = extractPos(results, steps)

pos = zeros(3*(steps+1),1);
j = 1;
k = 1;
for i = 1:steps+1
    pos(j:3*i) = results(k:k+2);
    j = j+3;
    k = k + 12;
end

pos = reshape(pos, 3, steps+1);
x_i = pos(1,:);
y_i = pos(2,:);
z_i = pos(3,:);
end