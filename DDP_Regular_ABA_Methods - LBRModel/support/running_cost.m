function l = running_cost(in1,u)
%RUNNING_COST
%    L = RUNNING_COST(IN1,U)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    05-Feb-2020 08:48:55

theta = in1(2,:);
xdot = in1(1,:);
l = theta.^2.*(1.1e1./8.0)+u.^2.*(3.0./1.0e1)+xdot.^2.*2.0;
