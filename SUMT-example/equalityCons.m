function h = equalityCons(U, x0, x_final_pos, info)
% equality constraints to be used in the ALM SUMT script
N = info.nMPC;
X = predictStates(x0, U, info);

h = (X((N)*12+1:N*12+3) - x_final_pos);

end