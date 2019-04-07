function x_k1 = forwardEuler(x,u,Ts,info)
%forwardEuler(x,u,Ts) takes the current state, input, and sampling time and
%computes the approximate discretization of states for use in an NLMPC 
%formulation.

f_xk = getStates(x,u,info);
x_k1 = x + f_xk*Ts;


end