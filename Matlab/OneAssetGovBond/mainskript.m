
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
casename            = 'SS_BASELINE_nm_100_fsolve';
% FindNewSS           = true;
FindNewSS           = false;
mpar.overrideEigen  = true;

%% Solve for Steady state
if FindNewSS

% Set parameters
defineSS_pars

% mainskript_steadystate
mainskript_SS_fsolve

else
    load(casename)
end
%%
SS_stats
%% Select aggregate shock
% aggrshock           = 'Uncertainty';
% par.rhoS            = 0.839;    % Persistence of variance
% par.sigmaS          = 0.539;    % STD of variance shocks
% aggrshock           = 'TFP';
% par.sigmaS          = 0.00965;
% par.rhoS            = 0.9;
aggrshock           = 'MP';
par.rhoS            = 0.0;    % Persistence of variance
par.sigmaS          = 0.001;    % STD of variance shocks

% par.prob_priceadj = 5/6;

% par.theta_pi    = 2.0; % Reaction to inflation
par.rho_R       = 0.0;  % Inertia

%% Produce matrices to reduce state-space
mpar.maxdim=100;

mainskript_statereduc

disp('Computing system for SS.');
toc
F = @(a,b,c,d)Fsys(a,b,c,d,Xss,Yss,Gamma_state,Gamma_control,InvGamma,Copula,par,mpar,grid,targets,P_H,aggrshock,oc);
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


plot_IRFs

