function [g] = inequalCon(U, Ts)
   % inequality contraints for the system
   % U(1) - Front, U(2) - Back, U(3) - Left, U(4) - Right
   F_upper = 10; T_lower = 0.1; T_upper = 500;
  
   g = [-U(1); U(1)/F_upper - 1;...
       -U(2); U(2)/F_upper - 1;...
       -U(3); U(3)/F_upper - 1;...
       -U(4); U(4)/F_upper - 1;...
       -Ts/T_lower + 1; Ts/T_upper - 1];
       
       
   