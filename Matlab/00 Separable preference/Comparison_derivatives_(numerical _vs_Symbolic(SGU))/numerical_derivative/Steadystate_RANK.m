

Gex = @(Guess)excessK_RANK(Guess,par);

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000000);

SS_return = fsolve(Gex,Guess,options);

grid.K = SS_return;

[excess,par,AggregateBondDemand,L,WW,Y]= excessK_RANK(SS_return,par);

grid.B = AggregateBondDemand;

% C=(par.RB-1)*grid.B + WW*L + par.PROFITS ;
C = Y - par.delta*grid.K;

% save RANK_lambda_35

function [excess,par,AggregateBondDemand,L,WW,Y]= excessK_RANK(Guess,par)

K = Guess;
% K = ( (par.RB+RkmR+par.delta-1)/(par.mu*(1-par.alpha)*(par.mu*par.alpha)^(par.alpha/(1-par.alpha+par.gamma))) )...
%        ^((1-par.alpha+par.gamma)/(-par.alpha*par.gamma));
L = (par.mu*par.alpha*K^(1-par.alpha))^(1/(1-par.alpha+par.gamma));
Y = L^(par.alpha)*K^(1-par.alpha);
par.Rk = par.mu*(1-par.alpha)*Y/K-par.delta+1;
WW = par.mu*par.alpha*Y/L;
par.L=L;
%% bank
RkmR=par.Rk-par.RB;
aa = par.lambda*par.beta*par.theta*RkmR;
bb = -par.lambda*(1-par.beta*par.theta*par.RB) + par.beta*(par.Rk-par.RB)*(1-par.theta) ;
cc = par.beta*par.RB*(1-par.theta);
  
par.phiB = ( -bb-sqrt(bb^2-4*aa*cc) )/(2*aa);
par.z = RkmR*par.phiB+par.RB;
par.x = RkmR*par.phiB+par.RB;
par.nuB = ( (1-par.theta)*par.beta*RkmR )/(1-par.beta*par.theta*par.x);
par.etaB = par.beta*par.RB*(1-par.theta)/(1-par.beta*par.theta*par.z); % beta*RB < 1 in incomlete market
  
par.NetWorth     =   par.omega*K/(1-par.theta*(RkmR*par.phiB+par.RB));
par.Ne    =   par.theta*par.z*par.NetWorth;
par.Nn    =   par.omega*K;

par.PROFITS = (1-par.mu)*Y;

AggregateCapitalDemand = par.phiB*par.NetWorth;
AggregateBondDemand = (par.phiB-1)*par.NetWorth;

excess=(K-AggregateCapitalDemand);
end
