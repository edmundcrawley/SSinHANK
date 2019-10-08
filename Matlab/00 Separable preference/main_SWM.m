%%% this solves ECON 5300, 2015, Ex 7.1(f)-(h)

%% 0. housekeeping
clear all;
clc;

%% 1. parameters and functional forms

% parameters

mmin   = 0;                  % borrowing constraint
par.beta   = 0.985;                % subjective discount factor
par.delta  = 0.02;                 % depreciation rate of physical capital
par.varphi = 2/3;                   % Frisch elasticity of labor supply
par.gamma  = 1/par.varphi;          % inverse elasticity of intertemporal substitution
par.eta=20;
par.mu=(par.eta-1)/par.eta;
par.alpha  = 2/3/par.mu;                   % labor income share
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)
par.PI=1;

rho    = 9.7/10;                  % persistence parameter prodictivity
alpha=1-par.alpha; varphi = par.varphi; delta = par.delta;
% set up asset grid

nh = 3;  mpar.nh=nh;              % number of possible productivity realizations
nm  = 100; mpar.nm=nm;            % number of asset grid points

mmax = 45;                        % maximum asset level

% A  = linspace(mmin,mmax,nm)';          % equally-spaced asset grid from a_1=b to a_M
grid.m = (exp(linspace(0,log(mmax - mmin+1),mpar.nm))-1+mmin);   
% grid.m = (exp(exp(exp(linspace(0,log(log(log(mmax - mmin+1)+1)+1),mpar.nm))-1)-1)-1+mmin);

A=grid.m';

    % Income Process
par.rhoH        = 0.979;    % Persistence of productivity
par.sigmaH      = 0.059;    % STD of productivity shocks
par.rhoS        = 0.839;    % Persistence of variance
par.sigmaS      = 0.539;    % STD of variance shocks

mpar.tauchen = 'importance';
mpar.in = 0.0005;
mpar.out = 0.0625;

%%%%%%%%%%%%%%%%%%%%%% grid h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P_H,grid,par]=stochastics_variance(par, mpar,grid);
[hgrid,P_H_,boundsH] = Tauchen(par.rhoH,mpar.nh,1, 0, mpar.tauchen); % LR variance = 1
% Correct long run variance for *human capital*
hgrid = hgrid*par.sigmaH/sqrt(1-par.rhoH^2);
hgrid = exp(hgrid); % Levels instead of Logs
grid.h=hgrid;
P_H = transition(mpar.nh,par.rhoH,sqrt(1-par.rhoH^2),boundsH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set up productivity grid
% H  = [h1,h2,h3]';                   % grid for productivity

H = grid.h;
% vectorize the grid in two dimensions
meshes.m = repmat(A,1,nh);            % values of A change vertically
meshes.h = repmat(H,nm,1);           % values of Y change horizontally

% % this is the built-in alternative
% [meshes.m,meshes.h] = ndgrid(A,Y);

% (inverse) marginal utility functions
up    = @(c) c.^(-par.gamma);        % marginal utility of consumption
invup = @(x) x.^(-1/par.gamma);      % inverse of marginal utility of consumption
vp    = @(h) h.^(1/par.varphi);      % marginal disutility of labor supply
invvp = @(x) x.^(par.varphi);        % inverse marginal disutility of labor supply    

%% 2. discretization



%% 3. endogenous functions

% optimal labor supply
L  = @(c,h,w) invvp(up(c).*h.*w);

% current consumption level, cp0(anext,ynext) is the guess
% C0 is equivalent to c_aux
C0 = @(cp0,r) invup(par.beta*(1+r)*up(cp0)*P_H');
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,h,c0,r,w) 1/(1+r)*(c0+anext-L(c0,h,w).*h.*w);
                
%% 4. solve for the stationary equilibrium

% convergence criterion for consumption iteration
crit = 10^(-10);


r0  = (1/par.beta-1);
% set up an anonymous function
fprintf('Start solving the Aiyagari model... \n');
tic;
myfun   = @(r) stationary_eqm_SWM(r0,crit,meshes,alpha,delta,rho,varphi,A0,C0,L,mpar,grid,P_H);
                                      
% options = optimset('display','iter','TolX',1e-8,'MaxIter',20);
options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000000);

rstar   = fsolve(myfun,r0,options);
% % [rstar,dist]   = fsolve(myfun,r0,options);

[cp0,c_update,L_new,jd]= stationary_eqm_SWM(r0,crit,meshes,alpha,delta,rho,varphi,A0,C0,L,mpar,grid,P_H);
                                                                
% w = par.mu*par.alpha*(par.mu*(1-par.alpha)/(rstar+par.delta))^((1-par.alpha)/par.alpha);
% 
% c0_ss = C0(cp0,rstar);
% 

% N_=L_new(:)'*jd(:);

% fprintf('Done with the Aiyagari model in %f sec. \n',toc);
% 
% % get the simulated asset levels
% fprintf('Fetching the wealth distribution... \n');
% [r,at] = stationary_equilibrium(rstar,crit,I,T,meshes.m,meshes.h,par.alpha,b,par.delta,rho,par.varphi,A0,C0,H);
% 
% %% 5. plot the wealth distribution
% 
% % use the last 100 periods
% [n,xout] = hist(at(:,T-100:T),M);     % choose M bins
% bar(xout,n/sum(n));                   % relative frequency is n/sum(n)
% title(sprintf('Wealth distribution, rho = %2.1f',rho));
% xlabel('Asset level')
% ylabel('Relative frequency')
