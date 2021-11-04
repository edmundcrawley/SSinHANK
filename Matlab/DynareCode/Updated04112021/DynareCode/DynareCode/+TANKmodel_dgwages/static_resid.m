function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = TANKmodel_dgwages.static_resid_tt(T, y, x, params);
end
residual = zeros(14, 1);
lhs = y(1);
rhs = params(2)*y(1)+T(1)*y(2);
residual(1) = lhs - rhs;
lhs = y(12);
rhs = params(4)*y(2)+params(5)*y(7);
residual(2) = lhs - rhs;
lhs = y(8);
rhs = y(9);
residual(3) = lhs - rhs;
lhs = y(10);
rhs = y(10)-1/params(4)*(y(6)-y(1));
residual(4) = lhs - rhs;
lhs = (1-params(11)*(1-params(2)))*y(11);
rhs = y(12)+y(9)-params(11)*(y(6)-y(1)-y(5))-y(5)*params(2)*params(11);
residual(5) = lhs - rhs;
lhs = y(2);
rhs = y(10)*params(13)+y(11)*params(12);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = y(8)*params(15)+y(9)*params(14);
residual(7) = lhs - rhs;
lhs = y(6);
rhs = y(1)*params(6)+y(2)*params(7)+y(13);
residual(8) = lhs - rhs;
lhs = y(5);
rhs = y(6)-y(1);
residual(9) = lhs - rhs;
residual(10) = y(3);
lhs = y(2);
rhs = y(4);
residual(11) = lhs - rhs;
lhs = y(13);
rhs = y(13)*params(3)+x(1);
residual(12) = lhs - rhs;
residual(13) = y(14);
lhs = y(4);
rhs = (1-params(1))*y(7);
residual(14) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
