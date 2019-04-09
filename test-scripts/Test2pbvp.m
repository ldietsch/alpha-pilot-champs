%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
addpath('../tpbvp');
g = 9.81;
Jx = 1; Jy = 1; Jz = 1;
m = 1; l_m = 1; k_tau = 1;
Mu1 = 0; Mu2 = 0; Mu3 = 0; Mu4 = 0; Mu5 = 0; Mu6 = 0; Mu7 = 0;

z = zeros(24,1);

solinit = bvpinit(linspace(0,15), z);
sol = bvp4c(@Quad_Model, @bc2, solinit);

x = sol.y(1,:);
y = sol.y(2,:);
z = sol.y(3,:);

t = sol.x;

L6 = sol.y(1,:);
L10 = sol.y(2,:);
L11 = sol.y(3,:);
L12 = sol.y(1,:);



figure
plot(t,x)
figure
plot(t,y)
figure
plot(t,z)

u1 = L6/m - Mu4 + (L12*k_tau)/Jz - (L11*l_m)/Jy;
u2 = L6/m - Mu5 - (L12*k_tau)/Jz + (L10*l_m)/Jx;
u3 = L6/m - Mu6 + (L12*k_tau)/Jz + (L11*l_m)/Jy;
u4 = L6/m - Mu7 - (L12*k_tau)/Jz - (L10*l_m)/Jx;

figure 
subplot(2,2,1)
plot(t,u1)

subplot(2,2,2)
plot(t,u2)

subplot(2,2,3)
plot(t,u3)

subplot(2,2,4)
plot(t,u4)

