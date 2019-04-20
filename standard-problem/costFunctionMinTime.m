function J = costFunctionMinTime(U, Q, R, x0, info)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the diagonal as Q0, N times1
rho = info.rho;
info.Ts = U(end);
X = predictStates(x0, U, info);
% X = X(13:end);
X = extractPosVec(X,info.nMPC);
X = X(4:end);
J = 1e-3*(X)'*Q*(X) + rho*U(1:end-1)'*R*U(1:end-1) + info.nMPC*...
    info.substeps*U(end);

end