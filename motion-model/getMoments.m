function J = getMoments(M,m, R, l)

J(1) = 2*M*R^2/5 + 2*l^2*m;
J(2) = 2*M*R^2/5 + 2*l^2*m;
J(3) = 2*M*R^2/5 + 4*l^2*m;

end