function [excess,cp0,joint_distr] = stationary_fsolve(r0,crit,Amat,Ymat,alpha,b,delta,rho,varphi,A0,C0,H,pi)
% this function
% (1) solves for the consumption decision rule, given an
% intereste rate, r0
% (2) simulates the stationary equilibrium associated with the interest
% rate, r0
% (3) (i) returns the residual between the given interes rate r0 and the one
% implied by the stationary aggregate capital and labor supply.
% (ii) returns as an optional output the wealth distribution.

% get dimensions of the grid
[M,N] = size(Amat);
mpar.nm=M;
mpar.nh=N;
% get productivity realizations from the first row
y1 = Ymat(1,1);
y2 = Ymat(1,2);

% compute the wage according to marginal pricing and a Cobb-Douglas 
% production function
w0 = (1-alpha)*(alpha/(r0+delta))^(alpha/(1-alpha));

% initial guess on future consumption (consume asset income plus labor
% income from working h=1.
cp0 = r0*Amat+Ymat*w0;

%%% iterate on the consumption decision rule

% distance between two successive functions
% start with a distance above the convergence criterion
dist    = crit+1;
% maximum iterations (avoids infinite loops)
maxiter = 10^(3);
% counter
iter    = 1;

fprintf('Inner loop, running... \n');
    
while (dist>crit&&iter<maxiter)

% derive current consumption
c0 = C0(cp0,r0);

% derive current assets
a0 = A0(Amat,Ymat,c0,r0,w0);

%%% update the guess for the consumption decision

% consumption decision rule for a binding borrowing constraint
% can be solved as a quadratic equation
% c^(2)-((1+r)a+b)c-(yw)^(1+2/3) = 0
% general formula and notation: ax^(2)+bx+c = 0
% x = (-b+sqrt(b^2-4ac))/2a
cpbind = ((1+r0)*Amat+b)/2+sqrt(((1+r0)*Amat+b).^(2)+4*(Ymat*w0).^(1+varphi))/2;
% % slow alternative: rootfinding gives the same result
% cpbind  = zeros(M,N);
% options = optimset('Display','off');
% cpbind(:,1) = fsolve(@(c) (1+r0)*Amat(:,1)+H(c,Y(1),w0).*Y(1).*w0-c+b,cp0(:,1),options);
% cpbind(:,2) = fsolve(@(c) (1+r0)*Amat(:,2)+H(c,Y(2),w0).*Y(2).*w0-c+b,cp0(:,2),options);

% consumption for nonbinding borrowing constraint
cpnon = zeros(M,N);
% interpolation conditional on productivity realization
% instead of extrapolation use the highest value in a0

for i=1:length(Ymat(1,:))
cpnon(:,i)  = interp1(a0(:,i),c0(:,i),Amat(:,i),'spline');
end

% merge the two, separate grid points that induce a binding borrowing constraint
% for the future asset level (the first observation of the endogenous current asset grid is the 
% threshold where the borrowing constraint starts binding, for all lower values it will also bind
% as the future asset holdings are monotonically increasing in the current
% asset holdings).

for i=1:length(Ymat(1,:))
cpnext(:,i) = (Amat(:,i)>a0(1,i)).*cpnon(:,i)+(Amat(:,i)<=a0(1,i)).*cpbind(:,i);
end

% distance measure
dist = norm((cpnext-cp0)./cp0);

% display every 100th iteration
% if mod(iter,100) == 1
%     fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dist]);
% end

% increase counter
iter = iter+1;

% update the guess on the consumption function
cp0 = cpnext;

end
fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dist]);
fprintf('Inner loop, done. \n');

m_star=a0;
P_H = pi;
grid.m = Amat(:,1)';
grid.h = Ymat(1,:);

[joint_distr]=JDiteration(m_star,P_H,mpar,grid);
  
joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);

AggregateSaving = m_star(:)'*joint_distr(:);
excess = 20 - AggregateSaving;

end

