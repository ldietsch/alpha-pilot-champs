function vEangles = getVelEangles(eAngles, bAngV)
phi = eAngles(1); theta = eAngles(2); %psi = eAngles(3);

vEangles = [1 sin(phi)*tan(theta) cos(phi)*tan(theta); ...
    0 cos(phi) -sin(phi);
    0 sin(phi)/cos(theta) cos(phi)/cos(theta)]*bAngV;


end