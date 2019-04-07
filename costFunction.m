function J = costFunction(x, xref, Q)
%The x's, xref's, and Q are for N steps ahead, so x is 12N x 1 and Q is PSD
%with the top left block a 3N x 3N identity matrix
J = (x - xref)'*Q*(x-xref);

end