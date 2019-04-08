function dxdt = Quad_Model(t, x)

g = 9.81;
Jx = 1; Jy = 1; Jz = 1;
m = 1; l_m = 1; k_tau = 1;
Mu1 = 1; Mu2 = 1; Mu3 = 1; Mu4 = 1; Mu5 = 1; Mu6 = 1; Mu7 = 1;

dxdt = [
x(6)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) - x(5)*(cos(x(7))*sin(x(9)) - cos(x(9))*sin(x(7))*sin(x(8))) + x(4)*cos(x(8))*cos(x(9))
x(5)*(cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(8))*sin(x(9))) - x(6)*(cos(x(9))*sin(x(7)) - cos(x(7))*sin(x(8))*sin(x(9))) + x(4)*cos(x(8))*sin(x(9))
x(4)*sin(x(8)) - x(6)*cos(x(7))*cos(x(8)) - x(5)*cos(x(8))*sin(x(7))
x(5)*x(12) - x(6)*x(11) - g*sin(x(8))
x(6)*x(10) - x(4)*x(12) + g*cos(x(8))*sin(x(7))
x(4)*x(11) - x(5)*x(10) + (Mu4 + Mu5 + Mu6 + Mu7 - (4*x(18))/m)/m + g*cos(x(7))*cos(x(8))
x(1) + x(12)*cos(x(7))*tan(x(8)) + x(11)*sin(x(7))*tan(x(8))
x(11)*cos(x(7)) - x(12)*sin(x(7))
(x(12)*cos(x(7)))/cos(x(8)) + (x(11)*sin(x(7)))/cos(x(8))
(x(11)*x(12)*(Jy - Jz))/Jx - (l_m*(Mu7 - Mu5 + (2*x(22)*l_m)/Jx))/Jx
- (l_m*(Mu4 - Mu6 + (2*x(23)*l_m)/Jy))/Jy - (x(10)*x(12)*(Jx - Jz))/Jy
(x(10)*x(11)*(Jx - Jy))/Jz - (k_tau*(Mu5 - Mu4 - Mu6 + Mu7 + (4*x(24)*k_tau)/Jz))/Jz
-x(19)
0
0
x(17)*x(12) - Mu1 - x(18)*x(11) - x(15)*sin(x(8)) - x(13)*cos(x(8))*cos(x(9)) - x(14)*cos(x(8))*sin(x(9))
x(18)*x(10) - x(16)*x(12) - Mu2 + x(13)*(cos(x(7))*sin(x(9)) - cos(x(9))*sin(x(7))*sin(x(8))) - x(14)*(cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(8))*sin(x(9))) + x(15)*cos(x(8))*sin(x(7))
x(16)*x(11) - Mu3 - x(17)*x(10) - x(13)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) + x(14)*(cos(x(9))*sin(x(7)) - cos(x(7))*sin(x(8))*sin(x(9))) + x(15)*cos(x(7))*cos(x(8))
x(20)*(x(12)*cos(x(7)) + x(11)*sin(x(7))) - x(21)*((x(11)*cos(x(7)))/cos(x(8)) - (x(12)*sin(x(7)))/cos(x(8))) - x(13)*(x(5)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) + x(6)*(cos(x(7))*sin(x(9)) - cos(x(9))*sin(x(7))*sin(x(8)))) + x(14)*(x(5)*(cos(x(9))*sin(x(7)) - cos(x(7))*sin(x(8))*sin(x(9))) + x(6)*(cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(8))*sin(x(9)))) + x(15)*(x(5)*cos(x(7))*cos(x(8)) - x(6)*cos(x(8))*sin(x(7))) - x(19)*(x(11)*cos(x(7))*tan(x(8)) - x(12)*sin(x(7))*tan(x(8))) - x(17)*g*cos(x(7))*cos(x(8)) + x(18)*g*cos(x(8))*sin(x(7))
x(16)*g*cos(x(8)) - x(19)*(x(11)*sin(x(7))*(tan(x(8))^2 + 1) + x(12)*cos(x(7))*(tan(x(8))^2 + 1)) - x(13)*(x(6)*cos(x(7))*cos(x(8))*cos(x(9)) - x(4)*cos(x(9))*sin(x(8)) + x(5)*cos(x(8))*cos(x(9))*sin(x(7))) - x(14)*(x(5)*cos(x(8))*sin(x(7))*sin(x(9)) - x(4)*sin(x(8))*sin(x(9)) + x(6)*cos(x(7))*cos(x(8))*sin(x(9))) - x(21)*((x(12)*cos(x(7))*sin(x(8)))/cos(x(8))^2 + (x(11)*sin(x(7))*sin(x(8)))/cos(x(8))^2) - x(15)*(x(4)*cos(x(8)) + x(6)*cos(x(7))*sin(x(8)) + x(5)*sin(x(7))*sin(x(8))) + x(18)*g*cos(x(7))*sin(x(8)) + x(17)*g*sin(x(7))*sin(x(8))
x(13)*(x(5)*(cos(x(7))*cos(x(9)) + sin(x(7))*sin(x(8))*sin(x(9))) - x(6)*(cos(x(9))*sin(x(7)) - cos(x(7))*sin(x(8))*sin(x(9))) + x(4)*cos(x(8))*sin(x(9))) - x(14)*(x(6)*(sin(x(7))*sin(x(9)) + cos(x(7))*cos(x(9))*sin(x(8))) - x(5)*(cos(x(7))*sin(x(9)) - cos(x(9))*sin(x(7))*sin(x(8))) + x(4)*cos(x(8))*cos(x(9)))
x(18)*x(5) - x(17)*x(6) + (x(23)*x(12)*(Jx - Jz))/Jy - (x(24)*x(11)*(Jx - Jy))/Jz
x(16)*x(6) - x(18)*x(4) - x(20)*cos(x(7)) - x(19)*sin(x(7))*tan(x(8)) - (x(21)*sin(x(7)))/cos(x(8)) - (x(22)*x(12)*(Jy - Jz))/Jx - (x(24)*x(10)*(Jx - Jy))/Jz
x(17)*x(4) - x(16)*x(5) + x(20)*sin(x(7)) - x(19)*cos(x(7))*tan(x(8)) - (x(21)*cos(x(7)))/cos(x(8)) - (x(22)*x(11)*(Jy - Jz))/Jx + (x(23)*x(10)*(Jx - Jz))/Jy
];

end