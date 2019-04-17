function J = costFunction(U, Uref, Xref, Q, R, x0, info)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the diagonal as Q0, N times1
rho = info.rho;
X = predictPosition(x0, U, info); %to save time, we will only predict position
X = X(4:end); %dont care about initial values in cost function, start at 4
J = (X - Xref)'*Q*(X - Xref) + rho*(U-Uref)'*R*(U-Uref);

end