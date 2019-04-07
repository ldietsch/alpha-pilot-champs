function xref = genRefTraj(x0,xf,intervals, N)
steps = N*intervals;
delX = (xf - x0)/steps;

xref = zeros(steps+1,3);
xref(1,:) = x0;
for i = 2:steps+1
    xref(i,:) = xref(i-1,:)+delX';
end

end