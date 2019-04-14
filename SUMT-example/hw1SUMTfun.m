function f = hw1SUMTfun(x)
%HW1SUMTFUN
%Objective function as part of the pseduo-objective function in homework 1
%part 2. a - alpha, the angle of attack in radians and b - the length of 
%the wingspan
b=x(2); a=x(3);
b = b/100;
f = b.^2.*(a.*6.0e1+pi).^2;

end
