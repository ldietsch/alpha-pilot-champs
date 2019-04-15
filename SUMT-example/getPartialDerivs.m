function del = getPartialDerivs(x0, U, f, info)
% getPartialDerivs finds the partial derivatives of any function to be used
% with the ALM SUMT method. Note x0 and info are unused explicitly but must
% be passed into the function since f is parameterized by x0 and info.

fplus = zeros(size(f(U),1),length(U));
fminus = zeros(size(f(U),1),length(U));
if size(f(U),1) == 1
    for i = 1:length(U)
        U(i) = U(i) + eps;
        fplus(i) = f(U);
        U(i) = U(i) - 2*eps;
        fminus(i) = f(U);
        U(i) = U(i) + 2*eps; %set back to normal for next round
    end
del = (fplus - fminus)./(2*eps);
else
    for i = 1:length(U)
        U(i) = U(i) + eps;
        fplus(:,i) = f(U);
        U(i) = U(i) - 2*eps;
        fminus(:,i) = f(U);
        U(i) = U(i) + 2*eps; %set back to normal for next round        
    end
del = (fplus - fminus)./(2*eps);
end

end
