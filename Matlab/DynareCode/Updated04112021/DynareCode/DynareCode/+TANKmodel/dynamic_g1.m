function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = TANKmodel.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(14, 20);
g1(1,4)=1;
g1(1,18)=(-params(2));
g1(1,5)=(-T(1));
g1(2,11)=(-params(5));
g1(2,13)=(-params(4));
g1(2,15)=1;
g1(3,12)=(-params(5));
g1(3,14)=(-params(4));
g1(3,15)=1;
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
g1(11,5)=1;
g1(11,7)=(-1);
g1(12,3)=(-params(3));
g1(12,16)=1;
g1(12,20)=(-1);
g1(13,17)=1;
g1(14,7)=1;
g1(14,10)=(-(1-params(1)));

end
