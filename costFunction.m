function J = costFunction(U, Xref, Q)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the diagonal as Q0, N times

J = (X - Xref)'*Q*(X - Xref);

end