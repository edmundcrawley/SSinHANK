%% 0. housekeeping
clear all;
clc;

%% 1. parameters and functional forms

% parameters
alpha  = 1/3;                   % capital income share
b      = 0;                     % borrowing constraint
beta   = 97/100;                % subjective discount factor
delta  = 5/100;                 % depreciation rate of physical capital
gamma  = 3/2;                   % inverse elasticity of intertemporal substitution
varphi = 2/3;                   % Frisch elasticity of labor supply
rho    = 6/10;                  % persistence parameter prodictivity
N_h    = 2;                     % number of possible productivity realizations
y1     = 95/100;                % low productivity realization
y2     = 105/100;               % high productivity realization

% transition probability matrix (productivity)
pi      = zeros(N_h,N_h);

% probabilities have to sum up to 1
pi(1,1) = rho;                  % transition from state 1 to state 1
pi(1,2) = 1-rho;                % transition from state 2 to state 1
pi(2,1) = 1-rho;                % transition from state 1 to state 2
pi(2,2) = rho;                  % transition from state 2 to state 2

% (inverse) marginal utility functions
up    = @(c) c.^(-gamma);        % marginal utility of consumption
invup = @(x) x.^(-1/gamma);      % inverse of marginal utility of consumption
vp    = @(h) h.^(1/varphi);      % marginal disutility of labor supply
invvp = @(x) x.^(varphi);        % inverse marginal disutility of labor supply    

%% 2. discretization

% set up asset grid
% M  = 100;                        % number of asset grid points
% aM = 45;                         % maximum asset level
% A  = linspace(b,aM,M)';          % equally-spaced asset grid from a_1=b to a_M

% set up productivity grid
% Y  = [y1,y2]';                   % grid for productivity

% % vectorize the grid in two dimensions
% Amat = repmat(A,1,N);            % values of A change vertically
% Ymat = repmat(Y',M,1);           % values of Y change horizontally

% Income Process

par.rhoH        = 0.979;    % Persistence of productivity
par.sigmaH      = 0.059;    % STD of productivity shocks
par.rhoS        = 0.839;    % Persistence of variance
par.sigmaS      = 0.539;    % STD of variance shocks

%
mpar.in         = 0.0005;  % Prob. to become entrepreneur
mpar.out        = 0.0625;   % Prob. to become worker again
mpar.tauchen    ='importance';

%
mpar.nm         = 100;
mpar.nh         = 4;
m_max = 45; m_min=0;
grid.m = (exp(exp(exp(linspace(0,log(log(log(m_max - m_min+1)+1)+1),mpar.nm))-1)-1)-1+m_min);

[P_H,grid,par]=stochastics_variance_wo_E(par, mpar,grid);

[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);
%% 3. endogenous functions

% optimal labor supply
H  = @(c,h,w) invvp(up(c).*h.*w);

% current consumption level, cp0(anext,ynext) is the guess
% C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');
C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*P_H);
                
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,h,c0,r,w) 1/(1+r)...
                    *(c0+anext-H(c0,h,w).*h.*w);

%% 4. solve for the stationary equilibrium

% convergence criterion for consumption iteration
crit = 10^(-6);

% parameters of the simulation

% r0  = (1/beta-1)-[10^(-12),10^(-4)];
r0  = (1/beta-1);
% set up an anonymous function
fprintf('Start solving the Aiyagari model... \n');
tic;
myfun   = @(r) stationary_fsolve(r,crit,meshes.m,meshes.h,alpha,b,delta,rho,varphi,A0,C0,H,P_H);
% options = optimset('display','iter','TolX',1e-8,'MaxIter',20);
options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000000);

rstar   = fsolve(myfun,r0,options);

% update with equilibrium rstar
[excess,cp0,jd]=stationary_fsolve(rstar,crit,meshes.m,meshes.h,alpha,b,delta,rho,varphi,A0,C0,H,P_H);

fprintf('Done with the Aiyagari model in %f sec. \n',toc);

w_ss = (1-alpha)*(alpha/(rstar+delta))^(alpha/(1-alpha));

c_ss = C0(cp0,rstar);
a_ss = A0(meshes.m,meshes.h,c_ss,rstar,w_ss);
ni_ss = H(c_ss,meshes.h,w_ss);

N = sum(sum(jd.*ni_ss));
K = sum(sum(jd.*a_ss));
Y = K^alpha*N^(1-alpha);

