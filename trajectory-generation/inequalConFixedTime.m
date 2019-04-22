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

vel = extractVel(xref,info.Nsteps);
vf = vel(end-9:end-7);

vmax = 54; % m/s
beta = 2; %m/s
if theta < pi/2
    vUpper = (norm(vf) - (1 - theta/(pi/2))*vmax - beta)/...
        ((1 - theta/(pi/2))*vmax - beta);
    vLower = -(norm(vf) + beta)/beta;
else
    vUpper = (norm(vf) - beta)/beta;
    vLower = -norm(vf);
end

omega = extractOmega(xref,info.Nsteps);
omf = omega(end-2:end);
omfU = (norm(omf) - pi)/pi;
omfL = -(norm(omf) - pi)/pi;

g = [gUpper; vUpper; omfU; gLower; vLower; omfL];

