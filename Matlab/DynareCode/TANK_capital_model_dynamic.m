function [residual, g1, g2, g3] = TANK_Capital_model_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(17, 1);
Omega__ = (1-params(1))/(1-params(1)+params(1)*params(8));
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+(1-params(1))*params(4));
lambda_gali__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda_gali__*(params(4)+(params(1)+params(5))/(1-params(1)));
lhs =y(5);
rhs =params(2)*y(22)+kappa__*y(6);
residual(1)= lhs-rhs;
lhs =y(16);
rhs =params(4)*y(14)+params(5)*y(12);
residual(2)= lhs-rhs;
lhs =y(16);
rhs =params(4)*y(15)+params(5)*y(13);
residual(3)= lhs-rhs;
lhs =y(14);
rhs =y(24)-1/params(4)*(y(10)-y(22));
residual(4)= lhs-rhs;
lhs =y(15)*(1-params(11)*(1-params(2)));
rhs =y(16)+y(13)-params(11)*(y(2)-y(5)-y(1))-params(2)*params(11)*y(9);
residual(5)= lhs-rhs;
lhs =y(8);
rhs =y(14)*params(13)+y(15)*params(12)+params(17)*y(20);
residual(6)= lhs-rhs;
lhs =y(11);
rhs =y(12)*params(15)+y(13)*params(14);
residual(7)= lhs-rhs;
lhs =y(20)*params(16);
rhs =y(19)-(1-params(16))*y(4);
residual(8)= lhs-rhs;
lhs =y(21);
rhs =params(18)*(y(19)-y(4));
residual(9)= lhs-rhs;
lhs =y(9)+y(21);
rhs =params(2)*(1-params(16))*y(26)+(1-params(2)*(1-params(16)))*(y(25)+y(23)-y(19));
residual(10)= lhs-rhs;
lhs =y(10);
rhs =y(5)*params(6)+y(6)*params(7)+y(17);
residual(11)= lhs-rhs;
lhs =y(9);
rhs =y(10)-y(22);
residual(12)= lhs-rhs;
lhs =y(7);
rhs =psi_n_ya__*(y(18)+params(1)*y(4));
residual(13)= lhs-rhs;
lhs =y(6);
rhs =y(8)-y(7);
residual(14)= lhs-rhs;
lhs =y(17);
rhs =params(3)*y(3)+x(it_, 1);
residual(15)= lhs-rhs;
residual(16) = y(18);
lhs =y(8);
rhs =params(1)*y(4)+y(18)+(1-params(1))*y(11);
residual(17)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(17, 27);

  %
  % Jacobian matrix
  %

  g1(1,5)=1;
  g1(1,22)=(-params(2));
  g1(1,6)=(-kappa__);
  g1(2,12)=(-params(5));
  g1(2,14)=(-params(4));
  g1(2,16)=1;
  g1(3,13)=(-params(5));
  g1(3,15)=(-params(4));
  g1(3,16)=1;
  g1(4,22)=(-(1/params(4)));
  g1(4,10)=1/params(4);
  g1(4,14)=1;
  g1(4,24)=(-1);
  g1(5,5)=(-params(11));
  g1(5,1)=(-params(11));
  g1(5,9)=params(2)*params(11);
  g1(5,2)=params(11);
  g1(5,13)=(-1);
  g1(5,15)=1-params(11)*(1-params(2));
  g1(5,16)=(-1);
  g1(6,8)=1;
  g1(6,14)=(-params(13));
  g1(6,15)=(-params(12));
  g1(6,20)=(-params(17));
  g1(7,11)=1;
  g1(7,12)=(-params(15));
  g1(7,13)=(-params(14));
  g1(8,4)=1-params(16);
  g1(8,19)=(-1);
  g1(8,20)=params(16);
  g1(9,4)=params(18);
  g1(9,19)=(-params(18));
  g1(9,21)=1;
  g1(10,9)=1;
  g1(10,23)=(-(1-params(2)*(1-params(16))));
  g1(10,25)=(-(1-params(2)*(1-params(16))));
  g1(10,19)=1-params(2)*(1-params(16));
  g1(10,21)=1;
  g1(10,26)=(-(params(2)*(1-params(16))));
  g1(11,5)=(-params(6));
  g1(11,6)=(-params(7));
  g1(11,10)=1;
  g1(11,17)=(-1);
  g1(12,22)=1;
  g1(12,9)=1;
  g1(12,10)=(-1);
  g1(13,7)=1;
  g1(13,18)=(-psi_n_ya__);
  g1(13,4)=(-(params(1)*psi_n_ya__));
  g1(14,6)=1;
  g1(14,7)=1;
  g1(14,8)=(-1);
  g1(15,3)=(-params(3));
  g1(15,17)=1;
  g1(15,27)=(-1);
  g1(16,18)=1;
  g1(17,8)=1;
  g1(17,11)=(-(1-params(1)));
  g1(17,18)=(-1);
  g1(17,4)=(-params(1));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],17,729);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],17,19683);
end
end
