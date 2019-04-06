function a_body = getAccBody(bAngV, v_body, eAngles, m, F, g)
p = bAngV(1); q = bAngV(2); r = bAngV(3);
u = v_body(1); v = v_body(2); w = v_body(3);

a_body = [r*v-q*w; p*w-r*u; q*u-p*v] + [-g*sin(theta); ...
    g*cos(theta)*sin(phi); g*cos(theta)*cos(phi)] + 1/m*[0;0;-F];

end