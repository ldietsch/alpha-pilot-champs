function g = hw1SUMTcon(x)
%HW1SUMTCON
%Function to compute the constraints g(1) and g(2) are derived from the
%lift coefficients; g(3) from the min. lift required; g(4) and g(5) from
%the length of the min/max chord at the root; g(6) and g(7) from the
%min/max of the span; g(8) and g(9) from the min/max of the angle-of-attack
% a - alpha, the angle of attack in radians, c_r - chord length of
%the root, and b - the length of the wingspan
%cr = cr/10; 
cr = x(1); b=x(2); a=x(3);
g = [(b.*pi.*(a.*6.0e1+pi).*(1.0./2.7e2))-cr;...
    (b.*pi.*(a.*6.0e1+pi).*(-1.0./2.1e2))+cr;...
    b.^2.*(a.*6.0e1+3.141592653589793).*-1.256943310501668e-3+1.0;...
    1.0-cr/0.8;cr/3 - 1; ...
    b.*(-1.0./8.0)+1.0; b.*(1.0./1.4e1)-1.0;...
    -(a.*3.6e1)./pi-1.0;(a.*1.8e1)./pi-1.0];
end
