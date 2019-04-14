function [f, del, H] = potential(u,E,A,L,theta)
%Function to calculate potential energy function, its gradient, and Hessian
%in the homework assignment for AAE 550
k1 = E*A(1)/L(1); k2 = E*A(2)/L(2); k3 = E*A(3)/L(3);
K = [(k1+k2)*cos(theta)^2 (k2-k1)*sin(theta)*cos(theta); ...
    (k2 - k1)*sin(theta)*cos(theta) ((k1+k2)*sin(theta)^2 + k3)];
p = [P*cos(90-theta); P*sin(90-theta)];
%f is PI(u), the potential energy function
f = (1/2)*[u(1) u(2)]*(K*[u(1);u(2)]) - p'*[u(1); u(2)];
%del is the gradient of PI(u)
del = [(k1+k2)*cos(theta)^2*u(1) + (k2-k1)*sin(theta)*cos(theta)*u(2)- ...
    p(1); ((k1+k2)*sin(theta)^2+k3)*u(2) + ...
    (k2 - k1)*sin(theta)*cos(theta)*u(1) - p(2)];
%H is the Hessian of PI(u)
H = [(k1+k2)*cos(theta)^2 (k2 - k1)*sin(theta)*cos(theta); ...
    (k2 - k1)*sin(theta)*cos(theta) ((k1+k2)*sin(theta)^2 + k3)];
    
end

