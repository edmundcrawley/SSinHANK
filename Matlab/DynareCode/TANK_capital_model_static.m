function [residual, g1, g2] = TANK_Capital_model_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 17, 1);

%
% Model equations
%

Omega__ = (1-params(1))/(1-params(1)+params(1)*params(8));
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+(1-params(1))*params(4));
lambda_gali__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda_gali__*(params(4)+(params(1)+params(5))/(1-params(1)));
lhs =y(1);
rhs =params(2)*y(1)+kappa__*y(2);
residual(1)= lhs-rhs;
lhs =y(12);
rhs =params(4)*y(10)+params(5)*y(8);
residual(2)= lhs-rhs;
lhs =y(12);
rhs =params(4)*y(11)+params(5)*y(9);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(10)-1/params(4)*(y(6)-y(1));
residual(4)= lhs-rhs;
lhs =y(11)*(1-params(11)*(1-params(2)));
rhs =y(12)+y(9)-params(11)*(y(6)-y(1)-y(5))-y(5)*params(2)*params(11);
residual(5)= lhs-rhs;
lhs =y(4);
rhs =y(10)*params(13)+y(11)*params(12)+params(17)*y(16);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =y(8)*params(15)+y(9)*params(14);
residual(7)= lhs-rhs;
lhs =y(16)*params(16);
rhs =y(15)-y(15)*(1-params(16));
residual(8)= lhs-rhs;
residual(9) = y(17);
lhs =y(5)+y(17);
rhs =y(17)*params(2)*(1-params(16))+(1-params(2)*(1-params(16)))*(y(12)+y(7)-y(15));
residual(10)= lhs-rhs;
lhs =y(6);
rhs =y(1)*params(6)+y(2)*params(7)+y(13);
residual(11)= lhs-rhs;
lhs =y(5);
rhs =y(6)-y(1);
residual(12)= lhs-rhs;
lhs =y(3);
rhs =psi_n_ya__*(y(14)+params(1)*y(15));
residual(13)= lhs-rhs;
lhs =y(2);
rhs =y(4)-y(3);
residual(14)= lhs-rhs;
lhs =y(13);
rhs =y(13)*params(3)+x(1);
residual(15)= lhs-rhs;
residual(16) = y(14);
lhs =y(4);
rhs =params(1)*y(15)+y(14)+(1-params(1))*y(7);
residual(17)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(17, 17);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(2);
  g1(1,2)=(-kappa__);
  g1(2,8)=(-params(5));
  g1(2,10)=(-params(4));
  g1(2,12)=1;
  g1(3,9)=(-params(5));
  g1(3,11)=(-params(4));
  g1(3,12)=1;
  g1(4,1)=(-(1/params(4)));
  g1(4,6)=1/params(4);
  g1(5,1)=(-params(11));
  g1(5,5)=(-(params(11)-params(2)*params(11)));
  g1(5,6)=params(11);
  g1(5,9)=(-1);
  g1(5,11)=1-params(11)*(1-params(2));
  g1(5,12)=(-1);
  g1(6,4)=1;
  g1(6,10)=(-params(13));
  g1(6,11)=(-params(12));
  g1(6,16)=(-params(17));
  g1(7,7)=1;
  g1(7,8)=(-params(15));
  g1(7,9)=(-params(14));
  g1(8,15)=(-(1-(1-params(16))));
  g1(8,16)=params(16);
  g1(9,17)=1;
  g1(10,5)=1;
  g1(10,7)=(-(1-params(2)*(1-params(16))));
  g1(10,12)=(-(1-params(2)*(1-params(16))));
  g1(10,15)=1-params(2)*(1-params(16));
  g1(10,17)=1-params(2)*(1-params(16));
  g1(11,1)=(-params(6));
  g1(11,2)=(-params(7));
  g1(11,6)=1;
  g1(11,13)=(-1);
  g1(12,1)=1;
  g1(12,5)=1;
  g1(12,6)=(-1);
  g1(13,3)=1;
  g1(13,14)=(-psi_n_ya__);
  g1(13,15)=(-(params(1)*psi_n_ya__));
  g1(14,2)=1;
  g1(14,3)=1;
  g1(14,4)=(-1);
  g1(15,13)=1-params(3);
  g1(16,14)=1;
  g1(17,4)=1;
  g1(17,7)=(-(1-params(1)));
  g1(17,14)=(-1);
  g1(17,15)=(-params(1));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],17,289);
end
end
