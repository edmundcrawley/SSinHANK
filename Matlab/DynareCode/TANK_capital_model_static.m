function [residual, g1, g2] = TANK_capital_model_static(y, x, params)
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

residual = zeros( 18, 1);

%
% Model equations
%

Omega__ = 1;
YoverC__ = 1.0/(params(12)+params(13));
InvoverC__ = params(17)/(params(12)+params(13));
psi_n_ya__ = (1+params(5))/(params(1)+params(5)+params(4)*YoverC__*(1-params(1)));
lambda_gali__ = (1-params(9))*(1-params(9)*params(2))/params(9)*Omega__;
kappa__ = lambda_gali__*(params(4)*YoverC__+(params(5)+params(1))/(1-params(1)));
lhs =y(1);
rhs =params(2)*y(1)+kappa__*y(2);
residual(1)= lhs-rhs;
lhs =y(13);
rhs =params(4)*(params(13)*y(10)+params(12)*y(11))/(params(12)+params(13))+params(5)*y(7);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =y(9);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(10)-1/params(4)*(y(6)-y(1));
residual(4)= lhs-rhs;
lhs =y(11)*(1-params(11)*(1-params(2)));
rhs =y(13)+y(9)-params(11)*(y(6)-y(1)-y(5))-y(5)*params(2)*params(11);
residual(5)= lhs-rhs;
lhs =y(4);
rhs =params(13)*y(10)+params(12)*y(11)+params(17)*y(17);
residual(6)= lhs-rhs;
lhs =y(7);
rhs =y(8)*params(15)+y(9)*params(14);
residual(7)= lhs-rhs;
lhs =y(12);
rhs =(params(13)*y(10)+params(12)*y(11))/(params(12)+params(13));
residual(8)= lhs-rhs;
lhs =y(17)*params(16);
rhs =y(16)-y(16)*(1-params(16));
residual(9)= lhs-rhs;
residual(10) = y(18);
lhs =y(5)+y(18);
rhs =y(18)*params(2)*(1-params(16))+(1-params(2)*(1-params(16)))*(y(13)+y(7)-y(16));
residual(11)= lhs-rhs;
lhs =y(6);
rhs =y(1)*params(6)+y(2)*params(7)+y(14);
residual(12)= lhs-rhs;
lhs =y(5);
rhs =y(6)-y(1);
residual(13)= lhs-rhs;
lhs =y(3);
rhs =psi_n_ya__*(y(15)+params(1)*y(16))+y(17)*params(4)*(1-params(1))*InvoverC__/(params(1)+params(5)+YoverC__*params(4)*(1-params(1)));
residual(14)= lhs-rhs;
lhs =y(2);
rhs =y(4)-y(3);
residual(15)= lhs-rhs;
lhs =y(14);
rhs =y(14)*params(3)+x(1);
residual(16)= lhs-rhs;
residual(17) = y(15);
lhs =y(4);
rhs =params(1)*y(16)+y(15)+(1-params(1))*y(7);
residual(18)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(18, 18);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(2);
  g1(1,2)=(-kappa__);
  g1(2,7)=(-params(5));
  g1(2,10)=(-(params(13)*params(4)/(params(12)+params(13))));
  g1(2,11)=(-(params(12)*params(4)/(params(12)+params(13))));
  g1(2,13)=1;
  g1(3,8)=1;
  g1(3,9)=(-1);
  g1(4,1)=(-(1/params(4)));
  g1(4,6)=1/params(4);
  g1(5,1)=(-params(11));
  g1(5,5)=(-(params(11)-params(2)*params(11)));
  g1(5,6)=params(11);
  g1(5,9)=(-1);
  g1(5,11)=1-params(11)*(1-params(2));
  g1(5,13)=(-1);
  g1(6,4)=1;
  g1(6,10)=(-params(13));
  g1(6,11)=(-params(12));
  g1(6,17)=(-params(17));
  g1(7,7)=1;
  g1(7,8)=(-params(15));
  g1(7,9)=(-params(14));
  g1(8,10)=(-(params(13)/(params(12)+params(13))));
  g1(8,11)=(-(params(12)/(params(12)+params(13))));
  g1(8,12)=1;
  g1(9,16)=(-(1-(1-params(16))));
  g1(9,17)=params(16);
  g1(10,18)=1;
  g1(11,5)=1;
  g1(11,7)=(-(1-params(2)*(1-params(16))));
  g1(11,13)=(-(1-params(2)*(1-params(16))));
  g1(11,16)=1-params(2)*(1-params(16));
  g1(11,18)=1-params(2)*(1-params(16));
  g1(12,1)=(-params(6));
  g1(12,2)=(-params(7));
  g1(12,6)=1;
  g1(12,14)=(-1);
  g1(13,1)=1;
  g1(13,5)=1;
  g1(13,6)=(-1);
  g1(14,3)=1;
  g1(14,15)=(-psi_n_ya__);
  g1(14,16)=(-(params(1)*psi_n_ya__));
  g1(14,17)=(-(params(4)*(1-params(1))*InvoverC__/(params(1)+params(5)+YoverC__*params(4)*(1-params(1)))));
  g1(15,2)=1;
  g1(15,3)=1;
  g1(15,4)=(-1);
  g1(16,14)=1-params(3);
  g1(17,15)=1;
  g1(18,4)=1;
  g1(18,7)=(-(1-params(1)));
  g1(18,15)=(-1);
  g1(18,16)=(-params(1));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],18,324);
end
end
