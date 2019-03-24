%% Initialize workspace and load directories
clear
clc
close all
Computername='NKIM';
starttime=clock;
set(0,'defaulttextinterpreter','latex')
addpath(genpath('functions'))
addpath(genpath('latex'))

%% Select options
% Name economy
casename            = 'Tax_rebate_B_over_Y_4';
% FindNewSS           = true;
FindNewSS           = false;
mpar.overrideEigen  = true;

%% Solve for Steady state
if FindNewSS

%% Parameters
% Household Parameters
par.beta        = 0.985;     % Discount factor
par.xi          = 2;          % CRRA
par.gamma       = 1;          % Inverse Frisch elasticity

% Income Process
par.rhoH        = 0.979;    % Persistence of productivity 
par.sigmaH      = 0.059;    % STD of productivity shocks
mpar.in         = 0.00046;  % Prob. to become entrepreneur
mpar.out        = 0.0625;   % Prob. to become worker again

% Firm Side Parameters
par.eta         = 20; par.mu = (par.eta-1)/par.eta;       % Markup 5%
par.alpha       = 1;%2/3/par.mu;  % Labor share 2/3

% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

% Central Bank Policy
par.theta_pi    = 2; % Reaction to inflation
par.rho_R       = 0.0;  % Inertia (removed here because it is set later)


% Tax Schedule
par.tau         = 0.75;   % Proportional tax on labor and profit income 

% Debt rule
par.rho_B       = 0.0;  % Autocorrelation (removed here because it is set later)
par.gamma_pi    = 1.5;   % 1.5, Reaction to inflation
par.gamma_T     = 0;%0.5075; % Reaction to tax revenue

par.BtoY = 8;
%% Returns
par.PI          = 1;          % Gross inflation
par.RB          = 1;          % Market clearing interest rate to be solved for
par.borrwedge   = 1.02^0.25-1; % Wedge on borrowing

%% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nh         = 4;
mpar.tauchen    ='importance';

grid.K = 1; 

%% Numerical Parameters
mpar.crit    = 1e-11;

%% Quadruble Log Grid
m_min = 0; %natural borrowing limit (adjust according to eq interest rate)
m_max = 150; % 150
grid.m = exp(exp(linspace(0,log(log(m_max - m_min+1)+1),mpar.nm))-1)-1+m_min;

%% Use Tauchen method to approximate state space
[P_H,grid,par]=stochastics_variance(par, mpar,grid);

%% Solve for steady state

[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);
[c_guess,m_star,joint_distr,W_fc,Profits_fc,Output,N,grid,excess,par,C_agg,C_ind,inc,X_agg]...
    =steadystate_fsolve(P_H,grid,mpar,par,meshes);

%%
SS_targets

%% Prepare state and controls
grid.B=sum(grid.m.*sum(joint_distr,2)');

% Calculate Marginal Values of Capital (k) and Liquid Assets(m)
RBRB = (par.RB+(meshes.m<0)*par.borrwedge)./par.PI;

% Liquid Asset
mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment
Vm = RBRB.*mutil_c; %take return on money into account
Vm = reshape(Vm,[mpar.nm mpar.nh]);

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
%% Calculate Sufficient Statistiscs
Suff_Stats

%% adjust policy parameter after finding a steady state
par.theta_pi    = 2.5;
par.rho_R       = 0;  % Inertia

par.rho_B       = 0.99;
par.gamma_pi    = 1.25;   % 1.5, Reaction to inflation
par.gamma_T     = 0;%0.5075; % Reaction to tax revenue



%SS_stats
%% Select aggregate shock
% aggrshock           = 'Uncertainty';
% par.rhoS            = 0.839;    % Persistence of variance
% par.sigmaS          = 0.539;    % STD of variance shocks
% aggrshock           = 'TFP';
% par.sigmaS          = 0.00965;
% par.rhoS            = 0.0;
aggrshock           = 'MP';
par.rhoS            = 0.0;    % Persistence of variance
par.sigmaS          = 0.001;    % STD of variance shocks

%% Produce matrices to reduce state-space
mpar.maxdim=50;

%% mainskript_statereduc
%% Initialize state and control vector
tic
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
    squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
    log(par.RB); 0];

% Yss=[invmutil(mutil_c(:));log(par.PI);log(Output);log(par.W);log(par.PROFITS);log(par.N);log(targets.B);log(C_agg);...
%      par.G;log(targets.T);log(X_agg);log(par.mu)];

Yss=[invmutil(mutil_c(:));log(par.PI);log(Output);log(par.W);log(par.PROFITS);log(par.N);log(targets.B);log(C_agg);...
     par.G;log(X_agg);log(par.mu)];
%%
Statereduc_w_Gov

%%
disp('Computing system for SS.');
toc
F = @(a,b,c,d)Fsys_Gov_new(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,Copula,par,mpar,grid,targets,P_H,aggrshock,oc);
tic
[Fss,LHS,RHS,Distr] = F(State,State_m,Contr,Contr_m);
toc

%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,oc,mpar,par);

%% Produce IRFs
x0=zeros(mpar.numstates,1);
x0(end)=par.sigmaS;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;
mpar.maxlag=40;

for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end
aux=(squeeze(sum(sum(joint_distr,2),1)));
scale.h=repmat([1; aux(end); aux(1:end-2)],[1 mpar.maxlag]);


plot_IRFs_w_Gov_new