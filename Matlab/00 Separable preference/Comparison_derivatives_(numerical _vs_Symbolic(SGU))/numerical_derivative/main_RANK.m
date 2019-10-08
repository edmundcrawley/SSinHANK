clear all;clc;
%%
tic
par.xi = 2;
par.beta=0.99;
par.gamma=1;
par.epsilon=20;
par.mu = (par.epsilon-1)/par.epsilon;
par.alpha = 0.67;
par.delta = 0.08/4;
par.phi=11.4;
% RkmR_init=0.0025;
par.RB=1/par.beta;
grid.K=56.8661;
Guess=grid.K;
par.tau=1;
par.rho_R = 0.8;
par.theta_pi=1.25;
mpar.overrideEigen=true;

par.lambda = 0.3815;
par.theta = 0.972;
par.omega = 0.002;
%% Compute SS values
Steadystate_RANK

par.Q=1;
par.PI=1;
Output=Y;
par.W=WW;
par.RBB=par.Rk;
grid.N=L;
%% Select aggregate shock
aggrshock           = 'TFP';
% aggrshock           = 'MP';
% aggrshock           = 'NW';

% par.sigmaS          = 0.01;
par.rhoS_TFP            = 0.95;

par.rhoS_MP            = 0.0;    % Persistence of variance
% par.sigmaS          = 0.001;    % STD of variance shocks

par.rhoS_NW            = 0.0;    % Persistence of variance
% par.sigmaS          = -0.01;    % STD of variance shocks

% par.phi=0.1;

% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

%%

par.CB = (1-par.theta)*par.z*par.NetWorth-par.Nn;
par.C = C;
par.VR = (C-grid.N^(1+par.gamma)/(1+par.gamma))^(-par.xi);
par.Lambda = 1;

Xss=[log(grid.K);log(par.Q);log(par.VR);log(par.NetWorth);log(par.phiB);log(par.RB);log(par.RB);0;0;0];

par.I = par.delta*grid.K;

Yss=[log(Output);log(grid.N);log(par.I);log(par.C);log(par.Ne);...
     log(par.Nn);log(par.mu);log(par.W);log(par.PROFITS);log(Output);...
     log(par.Lambda); log(par.RBB);log(par.nuB);log(par.etaB);log(par.z);...
     log(par.x);log(par.PI)];

 
 
n1 = 0; % used for controls
n2 = 0; %used for distributions

mpar.nm=0;mpar.nh=0;
% Produce matrices to reduce state-space
oc = length(Yss);
os = length(Xss);

InvGamma                            = sparse(os+oc,os+oc);

% Gamma_state                                   = sparse(Gamma);
InvGamma(1:os,1:os)                 = eye(os);

Gamma_control                                 = sparse(oc,oc);
% Gamma_control(1:n1(1),1:n1(2))                = Poly;
% InvGamma(n2(2)+os+(1:n1(1)),n2(2)+os+(1:n1(2))) = InvCheb';

Gamma_control((1:oc),(1:oc))        = eye(oc);
InvGamma(os+(1:oc),os+(1:oc))       = eye(oc);

InvGamma                            = InvGamma';

mpar.numstates   = os;
mpar.numcontrols = oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;

%%
% F = @(a,b,c,d)Fsys_RANK_SGU_order(a,b,c,d,Xss,Yss,Gamma_control,InvGamma,par,mpar,grid,aggrshock,oc,os);
F = @(a,b,c,d)Fsys_RANK_SGU_order(a,b,c,d,Xss,Yss,Gamma_control,InvGamma,par,mpar,grid,aggrshock,oc,os);
tic
[Fss,LHS,RHS] = F(State,State_m,Contr,Contr_m);
toc

%% Solve RE via Schmitt-Grohe-Uribe Form
tic
[hx,gx,F1,F2,F3,F4,par] = SGU_solver_RANK(F, mpar.numstates,mpar.numcontrols,oc,mpar,par,grid);
toc


% x0=zeros(mpar.numstates,1);
x0=zeros(10,1); eta0 = zeros(mpar.numstates,1);
switch(aggrshock)
  case('TFP')
x0(end-2) = 0.01; % TFP shock
eta0(end-2) = 1;
  case('MP')
x0(end-1) = 0.001; % MP shock
eta0(end-1) = 1;
  case('NW')
x0(end) = -0.01; % TFP shock
eta0(end) = 1;
end

Levintal=true;

if Levintal==1
   tic
   [gxx_gss,hxx_hss,gss_L,hss_L,F11,F12,F13,F14,F21,F22,F23,F24,F31,F32,F33,F34,F41,F42,F43,F44,par] = ...
       SGU_solver_2nd_L(F,F1,F2,F3,F4,gx,hx,mpar.numstates,mpar.numcontrols,oc,mpar,par,grid,eta0);
   toc

   
   for i=1:mpar.numstates
    gxx_L(:,:,i) = gxx_gss(:,(i-1)*mpar.numstates+i:i*mpar.numstates+i-1);
    hxx_L(:,:,i) = hxx_hss(:,(i-1)*mpar.numstates+i:i*mpar.numstates+i-1);
   end
   
  gss_L = gxx_gss(:,end);
  hss_L = hxx_hss(:,end);
  
else
   tic
   [gxx,hxx,gss,hss,F11,F12,F13,F14,F21,F22,F23,F24,F31,F32,F33,F34,F41,F42,F43,F44,par] = ...
       SGU_solver_2nd(F,F1,F2,F3,F4,gx,hx,mpar.numstates,mpar.numcontrols,oc,mpar,par,grid,eta0);
   toc
end

toc
% MX=[eye(length(x0));gx];
% IRF_state_sparse=[];
% x=x0;
% mpar.maxlag=40;
% 
% for t=1:mpar.maxlag
%     IRF_state_sparse(:,t)=(MX*x)';
%     x=hx*x;
% end
% 
% plot_IRFs_w_Q