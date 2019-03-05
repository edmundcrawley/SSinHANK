%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for TANK model in Dynare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('c:\dynare'))
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First do a RANK model
dynare 'RepAgent.mod' noclearall;
RANK_irfs = oo_.irfs;

% MPC in this representative agent model is (1-beta)/(1+sigma/phi*(1-alpha))
% Note actually want MPC_hat = MPC/(MPC+MPS) = (1-beta) ? check this
MPC_RANK = (1-beta);
% Size of two channels (no other channels in rep agent model)
% 1) Aggregate income elasticity
Inc_wt_MPC_RANK = MPC_RANK;
% 2) Intertemporal Substitution Channel
Hicks_scaling_RANK = 1.0-MPC_RANK;
Elas_EIS_RANK = Hicks_scaling_RANK*sigma;

% Get inputs for partial eq. decomposition
% dY/Y
dY_Y = RANK_irfs.y_gap_eps_nu(1);
% dR/R
dR_R = RANK_irfs.r_real_eps_nu(1);

% Get output for partial eq. decomposition
% dC/C
dC_C = RANK_irfs.y_gap_eps_nu(1);

% Check they add up...
dC_C_Auclert = Inc_wt_MPC_RANK*dY_Y - Elas_EIS_RANK*dR_R;
error_RANK = dC_C_Auclert - dC_C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now the basic TANK model
sigma = 1;
phi=1;
phi_pi = 1.5;
phi_y  = .5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.99;
alpha=0.0;%1/3;
epsilon=6;
lambda = 0.3;
Lambda = 0.5;
% Calc steady state share of labor and consumption of each type
cons_share_to_labor_share_K = (1-Lambda*(1-beta))*(epsilon-1)/epsilon*(1-alpha);
cons_share_K_obj = @(x)x^sigma * (x/cons_share_to_labor_share_K)^phi - (lambda/(1-lambda))^(sigma+phi)*(1-x)^sigma * (1-x/cons_share_to_labor_share_K)^phi;
cons_share_K = fsolve(cons_share_K_obj, lambda);
%cons_share_K =cons_share_K(1);
cons_share_R = 1-cons_share_K;
labor_share_K = cons_share_K/cons_share_to_labor_share_K;
labor_share_R = 1-labor_share_K;
dynare 'TANKmodel.mod' noclearall;
TANK_irfs = oo_.irfs;

% Now calc Auclert's statistics
MPC_TANK_K = 1.0;
MPC_TANK_R = 1-beta;
MPC_TANK = MPC_TANK_R*cons_share_R + MPC_TANK_K*cons_share_K;
% 1) Aggregate income elasticity
Inc_wt_MPC_TANK = MPC_TANK;
% 2) Income Heterogeneity Channel
% Calculate this exactly (rather than finding gamma)
y_K = TANK_irfs.w_real_eps_nu(1) + TANK_irfs.n_K_eps_nu(1);
y_R = TANK_irfs.w_real_eps_nu(1) + TANK_irfs.n_R_eps_nu(1);
IncomeChannel_Total = y_K*cons_share_K + y_R*(1-beta)*cons_share_R;
% 3) Interest Rate Exposure Channel (URE measures as a fraction of total
% consumption/income)
y_K_share = 1/(1-Lambda*(1-beta))*cons_share_K; %share of total income earned by Keynesian labor
URE_K = -beta*Lambda*y_K_share;  % Nominal debt of Keynesians is a multiple of their labor income
URE_R = - URE_K;
Elas_R_TANK = URE_K*MPC_TANK_K + URE_R*MPC_TANK_R;
% 4) Fisher channel (debt deflation)
NNP_K = -Lambda*y_K_share;
NNP_R = - NNP_K;
Elas_P_TANK = NNP_K*MPC_TANK_K + NNP_R*MPC_TANK_R;
% 5) Intertemporal Substitution Channel
Hicks_scaling_TANK = (1.0-MPC_TANK_R)*cons_share_R + (1.0-MPC_TANK_K)*cons_share_K;
Elas_EIS_TANK = Hicks_scaling_TANK*sigma;

% Get inputs for partial eq. decomposition
% dY/Y
dY_Y_TANK = TANK_irfs.y_gap_eps_nu(1);
dYK_Y_TANK = (TANK_irfs.w_real_eps_nu(1) + TANK_irfs.n_K_eps_nu(1))*y_K_share;
dYR_Y_TANK = dY_Y_TANK - dYK_Y_TANK;
% dR/R
dR_R_TANK = TANK_irfs.r_real_eps_nu(1);
% dP/P
dP_P_TANK = TANK_irfs.pi_eps_nu(1);

% Get output for partial eq. decomposition
% dC/C
dC_C_TANK = TANK_irfs.y_gap_eps_nu(1);

% Check they add up...
dC_C_Auclert_TANK = MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK ...
                    + Elas_R_TANK*dR_R_TANK - Elas_P_TANK*dP_P_TANK ... 
                    - Elas_EIS_TANK*dR_R_TANK;
error_TANK = dC_C_Auclert_TANK - dC_C_TANK;

%check for Keynesians
dC_K = MPC_TANK_K*dYK_Y_TANK + URE_K*MPC_TANK_K*dR_R_TANK ...
        - NNP_K*MPC_TANK_K*dP_P_TANK;
dC_K- TANK_irfs.c_K_eps_nu(1)*cons_share_K

%check for Ricardians
dC_R = MPC_TANK_R*dYR_Y_TANK + URE_R*MPC_TANK_R*dR_R_TANK ...
        - NNP_R*MPC_TANK_R*dP_P_TANK ...
        - Elas_EIS_TANK*dR_R_TANK;
dC_R- TANK_irfs.c_R_eps_nu(1)*(1-cons_share_K)

