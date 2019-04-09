function vi = getVelInert(eAngles, v_body)
phi = eAngles(1);
theta = eAngles(2);
psi = eAngles(3);

DCM = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) ...
    cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) ...
    cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi); ...
    sin(theta) -sin(phi)*cos(theta) -cos(phi)*cos(theta)];

vi = DCM*v_body;

end