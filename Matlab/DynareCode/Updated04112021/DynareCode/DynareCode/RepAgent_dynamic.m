function [residual, g1, g2, g3] = RepAgent_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
Omega__ = 1;
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+params(4)*(1-params(1)));
lambda__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda__*(params(4)+(params(5)+params(1))/(1-params(1)));
T39 = (-1)/params(4);
lhs =y(2);
rhs =params(2)*y(12)+kappa__*y(3);
residual(1)= lhs-rhs;
lhs =y(3);
rhs =T39*(y(8)-y(12)-y(6))+y(13);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =y(2)*params(6)+y(3)*params(7)+y(10);
residual(3)= lhs-rhs;
lhs =y(6);
rhs =params(4)*psi_n_ya__*(y(14)-y(11));
residual(4)= lhs-rhs;
lhs =y(7);
rhs =y(8)-y(12);
residual(5)= lhs-rhs;
lhs =y(4);
rhs =psi_n_ya__*y(11);
residual(6)= lhs-rhs;
lhs =y(3);
rhs =y(5)-y(4);
residual(7)= lhs-rhs;
lhs =y(10);
rhs =params(3)*y(1)+x(it_, 1);
residual(8)= lhs-rhs;
residual(9) = y(11);
lhs =y(5);
rhs =y(11)+(1-params(1))*y(9);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 15);

  %
  % Jacobian matrix
  %

  g1(1,2)=1;
  g1(1,12)=(-params(2));
  g1(1,3)=(-kappa__);
  g1(2,12)=T39;
  g1(2,3)=1;
  g1(2,13)=(-1);
  g1(2,6)=T39;
  g1(2,8)=(-T39);
  g1(3,2)=(-params(6));
  g1(3,3)=(-params(7));
  g1(3,8)=1;
  g1(3,10)=(-1);
  g1(4,6)=1;
  g1(4,11)=params(4)*psi_n_ya__;
  g1(4,14)=(-(params(4)*psi_n_ya__));
  g1(5,12)=1;
  g1(5,7)=1;
  g1(5,8)=(-1);
  g1(6,4)=1;
  g1(6,11)=(-psi_n_ya__);
  g1(7,3)=1;
  g1(7,4)=1;
  g1(7,5)=(-1);
  g1(8,1)=(-params(3));
  g1(8,10)=1;
  g1(8,15)=(-1);
  g1(9,11)=1;
  g1(10,5)=1;
  g1(10,9)=(-(1-params(1)));
  g1(10,11)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,225);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,3375);
end
end
