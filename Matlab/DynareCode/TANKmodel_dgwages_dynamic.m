function [residual, g1, g2, g3] = TANKmodel_dgwages_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(14, 1);
Omega__ = 1;
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+params(4)*(1-params(1)));
lambda_gali__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda_gali__*(params(4)+(params(5)+params(1))/(1-params(1)));
lhs =y(4);
rhs =params(2)*y(18)+kappa__*y(5);
residual(1)= lhs-rhs;
lhs =y(15);
rhs =params(4)*y(5)+params(5)*y(10);
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(12);
residual(3)= lhs-rhs;
lhs =y(13);
rhs =y(19)-1/params(4)*(y(9)-y(18));
residual(4)= lhs-rhs;
lhs =(1-params(11)*(1-params(2)))*y(14);
rhs =y(15)+y(12)-params(11)*(y(2)-y(4)-y(1))-params(2)*params(11)*y(8);
residual(5)= lhs-rhs;
lhs =y(5);
rhs =y(13)*params(13)+y(14)*params(12);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(11)*params(15)+y(12)*params(14);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(4)*params(6)+y(5)*params(7)+y(16);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =y(9)-y(18);
residual(9)= lhs-rhs;
lhs =y(6);
rhs =psi_n_ya__*y(17);
residual(10)= lhs-rhs;
lhs =y(5);
rhs =y(7)-y(6);
residual(11)= lhs-rhs;
lhs =y(16);
rhs =params(3)*y(3)+x(it_, 1);
residual(12)= lhs-rhs;
residual(13) = y(17);
lhs =y(7);
rhs =y(17)+(1-params(1))*y(10);
residual(14)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(14, 20);

  %
  % Jacobian matrix
  %

  g1(1,4)=1;
  g1(1,18)=(-params(2));
  g1(1,5)=(-kappa__);
  g1(2,5)=(-params(4));
  g1(2,10)=(-params(5));
  g1(2,15)=1;
  g1(3,11)=1;
  g1(3,12)=(-1);
  g1(4,18)=(-(1/params(4)));
  g1(4,9)=1/params(4);
  g1(4,13)=1;
  g1(4,19)=(-1);
  g1(5,4)=(-params(11));
  g1(5,1)=(-params(11));
  g1(5,8)=params(2)*params(11);
  g1(5,2)=params(11);
  g1(5,12)=(-1);
  g1(5,14)=1-params(11)*(1-params(2));
  g1(5,15)=(-1);
  g1(6,5)=1;
  g1(6,13)=(-params(13));
  g1(6,14)=(-params(12));
  g1(7,10)=1;
  g1(7,11)=(-params(15));
  g1(7,12)=(-params(14));
  g1(8,4)=(-params(6));
  g1(8,5)=(-params(7));
  g1(8,9)=1;
  g1(8,16)=(-1);
  g1(9,18)=1;
  g1(9,8)=1;
  g1(9,9)=(-1);
  g1(10,6)=1;
  g1(10,17)=(-psi_n_ya__);
  g1(11,5)=1;
  g1(11,6)=1;
  g1(11,7)=(-1);
  g1(12,3)=(-params(3));
  g1(12,16)=1;
  g1(12,20)=(-1);
  g1(13,17)=1;
  g1(14,7)=1;
  g1(14,10)=(-(1-params(1)));
  g1(14,17)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,400);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,8000);
end
end
