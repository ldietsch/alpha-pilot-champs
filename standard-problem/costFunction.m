function J = costFunction(U, Xref, Q, R, x0, info)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the diagonal as Q0, N times1
rho = info.rho;
X = predictStates(x0, U, info);
X = X(13:end);
J = (X - Xref)'*Q*(X - Xref) + rho*U'*R*U;

end