function [residual, g1, g2] = RepAgent_static(y, x, params)
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

residual = zeros( 10, 1);

%
% Model equations
%

Omega__ = 1;
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+params(4)*(1-params(1)));
lambda__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda__*(params(4)+(params(5)+params(1))/(1-params(1)));
T38 = (-1)/params(4);
lhs =y(1);
rhs =params(2)*y(1)+kappa__*y(2);
residual(1)= lhs-rhs;
lhs =y(2);
rhs =y(2)+T38*(y(7)-y(1)-y(5));
residual(2)= lhs-rhs;
lhs =y(7);
rhs =y(1)*params(6)+y(2)*params(7)+y(9);
residual(3)= lhs-rhs;
residual(4) = y(5);
lhs =y(6);
rhs =y(7)-y(1);
residual(5)= lhs-rhs;
lhs =y(3);
rhs =psi_n_ya__*y(10);
residual(6)= lhs-rhs;
lhs =y(2);
rhs =y(4)-y(3);
residual(7)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(3)+x(1);
residual(8)= lhs-rhs;
residual(9) = y(10);
lhs =y(4);
rhs =y(10)+(1-params(1))*y(8);
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(2);
  g1(1,2)=(-kappa__);
  g1(2,1)=T38;
  g1(2,5)=T38;
  g1(2,7)=(-T38);
  g1(3,1)=(-params(6));
  g1(3,2)=(-params(7));
  g1(3,7)=1;
  g1(3,9)=(-1);
  g1(4,5)=1;
  g1(5,1)=1;
  g1(5,6)=1;
  g1(5,7)=(-1);
  g1(6,3)=1;
  g1(6,10)=(-psi_n_ya__);
  g1(7,2)=1;
  g1(7,3)=1;
  g1(7,4)=(-1);
  g1(8,9)=1-params(3);
  g1(9,10)=1;
  g1(10,4)=1;
  g1(10,8)=(-(1-params(1)));
  g1(10,10)=(-1);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
end
end
