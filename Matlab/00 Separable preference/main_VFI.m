
%% 0. housekeeping
clear all;
clc;
Computername='NKIM';
starttime=clock;
set(0,'defaulttextinterpreter','latex')


%% Select options
% Name economy
% casename            = 'BtoY3_min0_max150_nm100_nh6_xi2_gam2_eqgrid';
% casename            = 'BtoY3_min0_max150_nm100_nh6_xi3_gam3_eqgrid';
casename            = 'BtoY3_min0_max150_nm100_nh6_xi10_gam10_eqgrid';
% FindNewSS           = true;
FindNewSS           = false;
mpar.overrideEigen  = true;

%% Solve for Steady state
if FindNewSS==1% parameters

mmin = 0;                  % borrowing constraint
par.beta   = 0.985;        % subjective discount factor
par.delta  = 0.02;         % depreciation rate of physical capital
par.varphi = 1/10;          % Frisch elasticity of labor supply
par.xi = 10;                % inverse elasticity of intertemporal substitution (RRA)
par.gamma  = 1/par.varphi;   
par.eta=20;
par.mu=(par.eta-1)/par.eta;
par.alpha  = 1;            % labor income share
par.prob_priceadj = 3/4;   % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)
par.PI=1;
par.borrwedge=1.00^0.25-1;
par.BtoY = 3;
par.tau=1;

% set up asset grid

nh = 4;  mpar.nh=nh;              % number of possible productivity realizations
nm  = 20; mpar.nm=nm;            % number of asset grid points

mmax = 150;                        % maximum asset level

% grid.m  = linspace(mmin,mmax,nm);          % equally-spaced asset grid from a_1=b to a_M
grid.m = (exp(linspace(0,log(mmax - mmin+1),mpar.nm))-1+mmin);   
% grid.m = exp(exp(linspace(0,log(log(mmax - mmin+1)+1),mpar.nm))-1)-1+mmin;
% grid.m = (exp(exp(exp(linspace(0,log(log(log(mmax - mmin+1)+1)+1),mpar.nm))-1)-1)-1+mmin);

A=grid.m';
grid.K=1;
% Income Process
par.rhoH        = 0.979;    % Persistence of productivity
par.sigmaH      = 0.059;    % STD of productivity shocks
par.rhoS        = 0.839;    % Persistence of variance
par.sigmaS      = 0.539;    % STD of variance shocks

mpar.tauchen = 'importance';
mpar.in = 0.00046;
mpar.out = 0.0625;

[P_H,grid,par]=stochastics_variance(par, mpar,grid);


% set up productivity grid
% H  = [h1,h2,h3]';                   % grid for productivity

H = grid.h;
% vectorize the grid in two dimensions
meshes.m = repmat(A,1,nh);            % values of A change vertically
meshes.h = repmat(H,nm,1);           % values of Y change horizontally
% meshes.h(:,1:3)=meshes.h(:,1:3)/par.H;
% % this is the built-in alternative
% [meshes.m,meshes.h] = ndgrid(A,Y);

% (inverse) marginal utility functions
up    = @(c) c.^(-par.xi);        % marginal utility of consumption
invup = @(x) x.^(-1/par.xi);      % inverse of marginal utility of consumption
vp    = @(h) h.^(par.gamma);      % marginal disutility of labor supply
invvp = @(x) x.^(1/par.gamma);    % inverse marginal disutility of labor supply    

%% 2. endogenous functions

% optimal labor supply
L  = @(c,h,w) invvp(up(c).*h.*w); % eq (3) in pdf

% current consumption level, cp0(anext,ynext) is the guess
% C0 is equivalent to c_aux
C0 = @(cp0,r) invup(par.beta*(1+r)*up(cp0)*P_H');
% current asset level, c0 = C0(cp0(anext,ynext))
A0 = @(anext,h,c0,r,w) 1/(1+r)*(c0+anext-L(c0,h,w).*h.*w);

% C0 and A0 are not used in the code, just for reference
%% 3. solve for the steady state

% convergence criterion for consumption iteration
crit = 10^(-11);


r0  = (1/par.beta-1);
% set up an anonymous function
fprintf('Start solving the Aiyagari model... \n');
tic;

par.RB=1.04^0.25; % initial guess for equilibrium interest rate

Guess = par.RB;

tic
[L_guess,jd_guess] = stationary_eqm_VFI(Guess,crit,meshes,par,L,mpar,grid,P_H);
toc

myfun   = @(r) stationary_eqm(r,L_guess,jd_guess,crit,meshes,par,L,mpar,grid,P_H);

options=optimset('Display','off','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',4000000);

Return_ss   = fsolve(myfun,Guess,options);


%% update using the equilibrium r

[excess,c_new,L_new,Y,joint_distr,m_star,inc] ...
    = stationary_eqm(Return_ss,L_guess,jd_guess,crit,meshes,par,L,mpar,grid,P_H);


par.mc   =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;

grid.Y=Y;
grid.B = m_star(:)'*joint_distr(:);
par.RB = Return_ss;
inc_int = grid.B*(par.RB-1);
grid.C=c_new(:)'*joint_distr(:);
% c_aux(:)'*joint_distr(:)
grid.N=L_new(:)'*joint_distr(:);
mutil_c = 1./(c_new.^par.xi); % marginal utility at consumption policy no adjustment
%% Produce non-parametric Copula
cum_dist = cumsum(cumsum(joint_distr,1),2);
marginal_m = cumsum(squeeze((sum(joint_distr,2))));
marginal_h = cumsum(squeeze((sum(joint_distr,1))));

Copula = griddedInterpolant({marginal_m,marginal_h},cum_dist,'spline');

%% Save
filename=[casename];
save(filename)    

else
    load(casename)
end

% Choose parameters of the policy rule
par.theta_pi=1.5;
par.rho_R = 0.0;

% par.phi = 11.4;

par.lambda_pi = 0;
par.spread = 0;

%%
solver = 'QZ';
% solver = 'Time_Itertation';
%% Select aggregate shock

% aggrshock      = 'TFP';
% par.sigmaS     = 0.01;    
% par.rhoS       = 0.99;
aggrshock           = 'MP';
par.rhoS            = 0.0;    % Persistence of variance
par.sigmaS          = 0.001;    % STD of variance shocks
% aggrshock           = 'NetWorth';
% par.rhoS            = 0.0;    % Persistence of variance
% par.sigmaS          = -0.01;    % STD of variance shocks

% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)


% Produce matrices to reduce state-space
mpar.maxdim=10;

%% mainskript_statereduc
% mainskript_statereduc
mainskript_statereduc_wo_Ni

%%
Suff_Stats
%%
% F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,Copula,par,mpar,grid,P_H,aggrshock,oc,os);
F = @(a,b,c,d)Fsys_wo_Ni(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,Copula,par,mpar,grid,P_H,aggrshock,oc,os);
                           
tic
[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);
toc
           
%% Solve RE via Schmitt-Grohe-Uribe Form
tic
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,oc,mpar,par); % F3:nfx, F1:nfxp, F4: nfy, F2:nfyp
toc

%% Produce IRFs
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;
mpar.maxlag=30;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end

aux=(squeeze(sum(sum(joint_distr,2),1)));
scale.h=repmat([1; aux(end); aux(1:end-2)],[1 mpar.maxlag]);

eta0=zeros(mpar.numstates,1);
eta0(end)=1;
x0_simul=zeros(mpar.numstates,1);

% plot_IRFs_one_fig_w_Ni
plot_IRFs_one_fig_wo_Ni
