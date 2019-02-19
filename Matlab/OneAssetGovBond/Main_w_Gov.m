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
casename            = 'SS_BASELINE_tau_05_w_Gov_gam1_xi1_B_over_4Y_50';
FindNewSS           = true;
% FindNewSS           = true;
mpar.overrideEigen  = true;

%% Solve for Steady state
if FindNewSS

%% Parameters
% Household Parameters
par.beta        = 0.985;     % Discount factor
par.xi          = 1;          % CRRA
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
par.theta_pi    = 1.25; % Reaction to inflation
% par.rho_R       = 0.0;  % Inertia (removed here because it is set later)


% Tax Schedule
par.tau         = 0.95;   % Proportional tax on labor and profit income 

% Debt rule
par.gamma_pi    = 1.5;   % 1.5, Reaction to inflation
par.gamma_T     = 0.5075; % Reaction to tax revenue
%par.rho_B       = 0.99;  % Autocorrelation (removed here because it is set later)

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

%%   MPCs
[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

% NW=par.gamma/(1+par.gamma).*(par.N/par.H).*par.W;
NW=(par.N/par.H).*par.W;
WW=NW*ones(mpar.nm,mpar.nh); %Wages
WW(:,end)=par.PROFITS*par.profitshare;

% MPC
WW_h=squeeze(WW(1,:)); WW_h_mesh=squeeze(WW(:,:).*meshes.h);
grid_h_aux=grid.h;
MPC_m = zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
%     MPC_m(:,hh)=gradient(squeeze(c_guess(:,hh)))./gradient(grid.m)';
    MPC_m(:,hh)=gradient(squeeze(C_ind(:,hh)))./gradient(grid.m)'; % MPC_m_ is same with MPC_m
end


MPC_h = zeros(mpar.nm,mpar.nh);

for mm=1:mpar.nm
    MPC_h(mm,:)=gradient(log(c_guess(mm,:)))./gradient(log(WW_h.*grid_h_aux));
%    MPC_h_(mm,:)=gradient(c_guess(mm,:))./gradient(WW_h.*grid_h_aux);
end

NNP = meshes.m;
URE = inc.labor + meshes.m - c_guess;
%URE = inc.labor + meshes.m - C_ind; %URE should use c_guess as inc.labor
%is the labor income including disutility of working, so need to take this
%away

MPC_m=MPC_m.*(WW_h_mesh./C_ind);    %???????? what is this for?
MPC_m=min(MPC_m,1); % prevent MPC_a larger than 1
%% Sufficient statistics in Auclert

% 1. Aggregate income channel: EI[Yi/Y MPC_m]

Inc_wt_MPC = sum(sum( WW_h_mesh/sum(sum(joint_distr.*WW_h_mesh)).*MPC_m.*joint_distr ));
% income includes interest rate income?

% 3. Fisher channel
Redist_elas_P = sum(sum(MPC_m.*NNP.*joint_distr)) - sum(sum(MPC_m.*joint_distr))*sum(sum(NNP.*joint_distr));

% 4. Interest rate exposure channel
Redist_elas_R = sum(sum(MPC_m.*URE.*joint_distr)) - sum(sum(MPC_m.*joint_distr))*sum(sum(URE.*joint_distr));

% 5. substitution channel
sig_i = par.xi^(-1)*c_guess./C_ind; % c_guess: composite goods consumption, C_ind: final goods consumption

Hick_scaling = sum(sum(sig_i.*(1-MPC_m).*C_ind.*joint_distr));


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
%%
par.rho_B       = 0.999;
par.rho_R       = 0.0;  % Inertia

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

Yss=[invmutil(mutil_c(:));log(par.PI);log(Output);log(par.W);log(par.PROFITS);log(par.N);log(targets.B);log(C_agg);...
     log(par.G);log(targets.T);log(X_agg);log(par.mu)];
%%
Statereduc_w_Gov

%%
disp('Computing system for SS.');
toc
F = @(a,b,c,d)Fsys_Gov(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,Copula,par,mpar,grid,targets,P_H,aggrshock,oc);
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


plot_IRFs_w_Gov