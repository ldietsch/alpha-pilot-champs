function g = inequalConFixedTime(U, xref, info)
Fupper = info.motorlimit;
gUpper = U./Fupper - 1;
gLower = -U;

g0 = info.startGate;
g1 = info.nextGate;
g2 = info.afterNextGate;
if ~isempty(g2)
    vec1 = g1 - g0;
    vec2 = g2 - g1;
    theta = acos(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
else
    theta = 0;
end

vf = xref(end-9:end-7);

vmax = 54; % m/s
beta = 0.01; %m/s
alpha = 0.4;
if theta < pi/2
    vUpper = (vf - (1/sqrt(3))*alpha*(1 - theta/(pi/2))*vmax-beta);
    vLower = (-vf+1/sqrt(3)*beta);
else
    vUpper = vf;
    vLower = -vf;
end

omega = extractOmega(xref,info.Nsteps+1);
omf = omega(end-2:end);
omfU = (norm(omf) - pi)/pi;
omfL = -(norm(omf) + pi)/pi;

g = [gUpper; vUpper; omfU; gLower; vLower; omfL];
% g = [gUpper; vUpper; gLower; vLower];
end
