function a_body = getAccBody(bAngV, v_body, eAngles, m, F, g)
p = bAngV(1); q = bAngV(2); r = bAngV(3);
u = v_body(1); v = v_body(2); w = v_body(3);
phi = eAngles(1); theta = eAngles(2);
if norm(v_body) > 0
    vEb = v_body/norm(v_body);
else
    vEb = [1;0;0];
end
b = 0.01;
drag = -b.*dot(v_body,v_body).*vEb;
a_body = [r*v-q*w; p*w-r*u; q*u-p*v] + [-g*sin(theta); ...
    g*cos(theta)*sin(phi); g*cos(theta)*cos(phi)] + 1/m*[0;0;-F] + drag;

end