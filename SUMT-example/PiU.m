function [f, del, H, p] = PiU(u)
%Function to calculate potential energy function, its gradient, and Hessian
%in the homework assignment for AAE 550
%Parameters to be used 
theta = asin(2.75*0.5/3); E = 1.73e6; %psi
A =[(0.65/2)^2*pi;(0.65/2)^2*pi;(0.8/2)^2*pi]; %in^2
L = [3*12;3*12; (2.75/2+(5-2.75))*12]; %in
P = 14000; %lbs

k1 = E*A(1)/L(1); k2 = E*A(2)/L(2); k3 = E*A(3)/L(3);
K = [(k1+k2)*cos(theta)^2 (k2-k1)*sin(theta)*cos(theta); ...
    (k2 - k1)*sin(theta)*cos(theta) ((k1+k2)*sin(theta)^2 + k3)];
p = [P*cos(58); P*sin(58)];
%f is PI(u), the potential energy function
f = (1/2)*[u(1) u(2)]*(K*[u(1);u(2)]) - p'*[u(1); u(2)];
% syms x1 x2
% fsyms = (1/2)*[x1 x2]*(K*[x1;x2]) - p'*[x1; x2];
% format long
% expand(fsyms)
%del is the gradient of PI(u)
del = [(k1+k2)*cos(theta)^2*u(1) + (k2-k1)*sin(theta)*cos(theta)*u(2)- ...
    p(1); ((k1+k2)*sin(theta)^2+k3)*u(2) + ...
    (k2 - k1)*sin(theta)*cos(theta)*u(1) - p(2)];
%H is the Hessian of PI(u)
H = [(k1+k2)*cos(theta)^2 (k2 - k1)*sin(theta)*cos(theta); ...
    (k2 - k1)*sin(theta)*cos(theta) ((k1+k2)*sin(theta)^2 + k3)];
    
end

