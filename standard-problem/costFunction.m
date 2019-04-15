function J = costFunction(U, Q, R, x0, info)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the diagonal as Q0, N times1
rho = info.rho;
X = predictStates(x0, U, info);
X = X(13:end);
J = (X)'*Q*(X) + rho*U(1:end-1)'*R*U(1:end-1)+U(end) + info.Nsteps*...
    info.substeps*U(end);

end