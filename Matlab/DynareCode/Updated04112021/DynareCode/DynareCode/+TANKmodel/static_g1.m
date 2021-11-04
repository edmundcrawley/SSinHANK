function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = TANKmodel.static_g1_tt(T, y, x, params);
end
g1 = zeros(14, 14);
g1(1,1)=1-params(2);
g1(1,2)=(-T(1));
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
g1(6,2)=1;
g1(6,10)=(-params(13));
g1(6,11)=(-params(12));
g1(7,7)=1;
g1(7,8)=(-params(15));
g1(7,9)=(-params(14));
g1(8,1)=(-params(6));
g1(8,2)=(-params(7));
g1(8,6)=1;
g1(8,13)=(-1);
g1(9,1)=1;
g1(9,5)=1;
g1(9,6)=(-1);
g1(10,3)=1;
g1(11,2)=1;
g1(11,4)=(-1);
g1(12,13)=1-params(3);
g1(13,14)=1;
g1(14,4)=1;
g1(14,7)=(-(1-params(1)));
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
