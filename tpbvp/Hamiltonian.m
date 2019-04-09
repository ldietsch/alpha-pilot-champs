%%%%  Hamiltonian %%%%%
clear; clc; close all

x = sym('X',[1 12]);
l = sym('L',[1 12]);
mu = sym('Mu',[1 7]);
u = sym('U',[1 4]);
ua = sym('U_a',[1 4]);
va = sym('V_a',[1 3]);
Jx = sym('Jx'); Jy = sym('Jy'); Jz = sym('Jz');
k = sym('k_tau'); L = sym('l_m'); g = sym('g');
m = sym('m');

Fx = [...
    cos(x(8))*cos(x(9))*x(4) + (sin(x(7))*sin(x(8))*cos(x(9)) - cos(x(7))*sin(x(9)))*x(5) + (cos(x(7))*sin(x(8))*cos(x(9)) + sin(x(7))*sin(x(9)))*x(6)
    cos(x(8))*sin(x(9))*x(4) + (sin(x(7))*sin(x(8))*sin(x(9)) + cos(x(7))*cos(x(9)))*x(5) + (cos(x(7))*sin(x(8))*sin(x(9)) - sin(x(7))*cos(x(9)))*x(6)
    sin(x(8))*x(4) - sin(x(7))*cos(x(8))*x(5) - cos(x(7))*cos(x(8))*x(6)
    (x(5)*x(12) - x(6)*x(11)) - g*sin(x(8))
    (x(6)*x(10) - x(4)*x(12)) + g*cos(x(8))*sin(x(7))
    (x(4)*x(11) - x(5)*x(10)) + g*cos(x(8))*cos(x(7)) - (u(1)+u(2)+u(3)+u(4))/m
    x(1)+sin(x(7))*tan(x(8))*x(11)+cos(x(7))*tan(x(8))*x(12)
    cos(x(7))*x(11) - sin(x(7))*x(12)
    (sin(x(7))/cos(x(8)))*x(11) + (cos(x(7))/cos(x(8)))*x(12)
    ((Jy-Jz)/Jx)*x(11)*x(12) + (L/Jx)*(u(4)-u(2))
    ((Jz-Jx)/Jy)*x(10)*x(12) + (L/Jy)*(u(1)-u(3))
    ((Jx-Jy)/Jz)*x(10)*x(11) + (k/Jz)*(u(2)+u(4)-u(1)-u(3))
    ];

C_eff = [x(4)-va(1) x(5)-va(2) x(6)-va(3) u(1)-ua(1) u(2)-ua(2) u(3)-ua(3) u(4)-ua(4)]';


H = 1 + l*Fx + mu*C_eff + .5*u*eye(4)*u.';


for i = 1:12
    lambda_dot(i,1) = -diff(H, x(i));
end

U_l = [];
for i = 1:4
    lambda(i,1) = diff(H, u(i));
    u(i) = solve(lambda(i,1),u(i));
end

Fx = [...
    cos(x(8))*cos(x(9))*x(4) + (sin(x(7))*sin(x(8))*cos(x(9)) - cos(x(7))*sin(x(9)))*x(5) + (cos(x(7))*sin(x(8))*cos(x(9)) + sin(x(7))*sin(x(9)))*x(6)
    cos(x(8))*sin(x(9))*x(4) + (sin(x(7))*sin(x(8))*sin(x(9)) + cos(x(7))*cos(x(9)))*x(5) + (cos(x(7))*sin(x(8))*sin(x(9)) - sin(x(7))*cos(x(9)))*x(6)
    sin(x(8))*x(4) - sin(x(7))*cos(x(8))*x(5) - cos(x(7))*cos(x(8))*x(6)
    (x(5)*x(12) - x(6)*x(11)) - g*sin(x(8))
    (x(6)*x(10) - x(4)*x(12)) + g*cos(x(8))*sin(x(7))
    (x(4)*x(11) - x(5)*x(10)) + g*cos(x(8))*cos(x(7)) - (u(1)+u(2)+u(3)+u(4))/m
    x(1)+sin(x(7))*tan(x(8))*x(11)+cos(x(7))*tan(x(8))*x(12)
    cos(x(7))*x(11) - sin(x(7))*x(12)
    (sin(x(7))/cos(x(8)))*x(11) + (cos(x(7))/cos(x(8)))*x(12)
    ((Jy-Jz)/Jx)*x(11)*x(12) + (L/Jx)*(u(4)-u(2))
    ((Jz-Jx)/Jy)*x(10)*x(12) + (L/Jy)*(u(1)-u(3))
    ((Jx-Jy)/Jz)*x(10)*x(11) + (k/Jz)*(u(2)+u(4)-u(1)-u(3))
    ];


