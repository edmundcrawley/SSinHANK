function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = TANKmodel_dgwages.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(14, 1);
lhs = y(4);
rhs = params(2)*y(18)+y(5)*T(1);
residual(1) = lhs - rhs;
lhs = y(15);
rhs = params(4)*y(5)+params(5)*y(10);
residual(2) = lhs - rhs;
lhs = y(11);
rhs = y(12);
residual(3) = lhs - rhs;
lhs = y(13);
rhs = y(19)-1/params(4)*(y(9)-y(18));
residual(4) = lhs - rhs;
lhs = (1-params(11)*(1-params(2)))*y(14);
rhs = y(15)+y(12)-params(11)*(y(2)-y(4)-y(1))-params(2)*params(11)*y(8);
residual(5) = lhs - rhs;
lhs = y(5);
rhs = y(13)*params(13)+y(14)*params(12);
residual(6) = lhs - rhs;
lhs = y(10);
rhs = y(11)*params(15)+y(12)*params(14);
residual(7) = lhs - rhs;
lhs = y(9);
rhs = y(4)*params(6)+y(5)*params(7)+y(16);
residual(8) = lhs - rhs;
lhs = y(8);
rhs = y(9)-y(18);
residual(9) = lhs - rhs;
residual(10) = y(6);
lhs = y(5);
rhs = y(7);
residual(11) = lhs - rhs;
lhs = y(16);
rhs = params(3)*y(3)+x(it_, 1);
residual(12) = lhs - rhs;
residual(13) = y(17);
lhs = y(7);
rhs = (1-params(1))*y(10);
residual(14) = lhs - rhs;

end
