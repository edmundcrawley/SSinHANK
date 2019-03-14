%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for TANK model in Dynare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
addpath(genpath('c:\dynare'))
clear all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%
% Now the basic TANK model
sigma = 1.0;
phi=1.0;
phi_pi = 1.5;
phi_y  = 0.0;%.5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.99;
alpha=0.33;
epsilon=6;
lambda = 0.3;
Lambda = 2.0;
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
Elas_EIS_TANK = Hicks_scaling_TANK/sigma;

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

% Agg income channel
MPC_TANK*dY_Y_TANK
% Heterogeneous Income Channel
MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK - MPC_TANK*dY_Y_TANK
% Unhedged Interest Rate Exposure
Elas_R_TANK*dR_R_TANK
% Fisher Channel
- Elas_P_TANK*dP_P_TANK
% Intertemporal Elasticity Channel
- Elas_EIS_TANK*dR_R_TANK

%check for Keynesians
dC_K = MPC_TANK_K*dYK_Y_TANK + URE_K*MPC_TANK_K*dR_R_TANK ...
        - NNP_K*MPC_TANK_K*dP_P_TANK;
dC_K- TANK_irfs.c_K_eps_nu(1)*cons_share_K

%check for Ricardians
dC_R = MPC_TANK_R*dYR_Y_TANK + URE_R*MPC_TANK_R*dR_R_TANK ...
        - NNP_R*MPC_TANK_R*dP_P_TANK ...
        - Elas_EIS_TANK*dR_R_TANK;
dC_R- TANK_irfs.c_R_eps_nu(1)*(1-cons_share_K)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Now the TANK_capital model
sigma = 1;
phi=1;
phi_pi = 1.5;
phi_y  = 0.0;%.5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.99;
alpha=0.33;
epsilon=6;
lambda = 0.3;
Lambda = 2.0;
delta = 0.025;
psi_c = 3;
% Calc steady state share of labor and consumption of each type
cons_share_to_labor_share_K = (1-Lambda*(1-beta))*(epsilon-1)/epsilon*(1-alpha);
invest_share = delta*alpha*(epsilon-1)/(epsilon*(1/beta - (1-delta)));
cons_share_K_obj = @(x)x^sigma * (x/cons_share_to_labor_share_K)^phi - (lambda/(1-lambda))^(sigma+phi)*(1-invest_share-x)^sigma * (1-x/cons_share_to_labor_share_K)^phi;
cons_share_K = fsolve(cons_share_K_obj, lambda);
cons_share_R = 1-invest_share-cons_share_K;
labor_share_K = cons_share_K/cons_share_to_labor_share_K;
labor_share_R = 1-labor_share_K;
dynare 'TANK_capital_model.mod' noclearall;
TANK_Capital_irfs = oo_.irfs;

% Now calc Auclert's statistics
MPC_TANK_Capital_K = 1.0;
MPC_TANK_Capital_R = 1-beta;
MPC_TANK_Capital = MPC_TANK_Capital_R*cons_share_R + MPC_TANK_Capital_K*(1-cons_share_K);
% 1) Aggregate income elasticity
Inc_wt_MPC_TANK_Capital = MPC_TANK_Capital;
% 2) Income Heterogeneity Channel
% ????
% 3) Interest Rate Exposure Channel (URE measures as a fraction of total
% consumption/income)
cons_share = cons_share_R+cons_share_K;
y_K_share = 1/(1-Lambda*(1-beta))*cons_share_K; %share of total income earned by Keynesian labor
URE_K = -beta*Lambda*y_K_share;  % Nominal debt of Keynesians is a multiple of their labor income
URE_R = - URE_K + invest_share;
Elas_R_TANK_Capital = (URE_K*MPC_TANK_Capital_K + URE_R*MPC_TANK_Capital_R)/cons_share;
% 4) Fisher channel (debt deflation)
NNP_K = -Lambda*y_K_share;
NNP_R = - NNP_K;
Elas_P_TANK_Capital = (NNP_K*MPC_TANK_Capital_K + NNP_R*MPC_TANK_Capital_R)/cons_share;
% 5) Intertemporal Substitution Channel
Hicks_scaling_TANK_Capital = (1.0-MPC_TANK_Capital_R)*cons_share_R/cons_share + (1.0-MPC_TANK_Capital_K)*cons_share_K/cons_share;
Elas_EIS_TANK_Capital = Hicks_scaling_TANK_Capital/sigma;

% Get inputs for partial eq. decomposition
% dY/Y
dY_Y_TANK_Capital = TANK_Capital_irfs.y_eps_nu(1);
dYK_Y_TANK_Capital = (TANK_Capital_irfs.w_real_eps_nu(1) + TANK_Capital_irfs.n_K_eps_nu(1))*y_K_share;
dYR_Y_TANK_Capital = dY_Y_TANK_Capital - dYK_Y_TANK_Capital;
% dR/R
dR_R_TANK_Capital = TANK_Capital_irfs.r_real_eps_nu(1);
% dP/P
dP_P_TANK_Capital = TANK_Capital_irfs.pi_eps_nu(1);

% Get output for partial eq. decomposition
% dC/C
dC_C_TANK_Capital = TANK_Capital_irfs.c_R_eps_nu(1)*cons_share_R/cons_share + TANK_Capital_irfs.c_K_eps_nu(1)*cons_share_K/cons_share;

% Check they add up...
dC_C_Auclert_TANK_Capital = (MPC_TANK_Capital_R*dYR_Y_TANK_Capital + MPC_TANK_Capital_K*dYK_Y_TANK_Capital)/cons_share ...
                    + Elas_R_TANK_Capital*dR_R_TANK_Capital - Elas_P_TANK_Capital*dP_P_TANK_Capital ... 
                    - Elas_EIS_TANK_Capital*dR_R_TANK_Capital;
error_TANK_Capital = dC_C_Auclert_TANK_Capital - dC_C_TANK_Capital;

% Agg income channel
(MPC_TANK_Capital*dY_Y_TANK_Capital)/cons_share
% Heterogeneous Income Channel
(MPC_TANK_Capital_R*dYR_Y_TANK_Capital + MPC_TANK_Capital_K*dYK_Y_TANK_Capital - MPC_TANK_Capital*dY_Y_TANK_Capital)/cons_share
% Unhedged Interest Rate Exposure
Elas_R_TANK_Capital*dR_R_TANK_Capital
% Fisher Channel
- Elas_P_TANK_Capital*dP_P_TANK_Capital
% Intertemporal Elasticity Channel
- Elas_EIS_TANK_Capital*dR_R_TANK_Capital
%Total
dC_C_Auclert_TANK_Capital

%check for Keynesians
dC_K = MPC_TANK_Capital_K*dYK_Y_TANK_Capital/cons_share ...
        + URE_K*MPC_TANK_Capital_K*dR_R_TANK_Capital/cons_share ...
        - NNP_K*MPC_TANK_Capital_K*dP_P_TANK_Capital/cons_share;
dC_K- TANK_Capital_irfs.c_K_eps_nu(1)*cons_share_K/cons_share

%check for Ricardians
dC_R = MPC_TANK_Capital_R*dYR_Y_TANK_Capital/cons_share ...
        + URE_R*MPC_TANK_Capital_R*dR_R_TANK_Capital/cons_share ...
        - NNP_R*MPC_TANK_Capital_R*dP_P_TANK_Capital/cons_share ...
        - Elas_EIS_TANK_Capital*dR_R_TANK_Capital;
dC_R- TANK_Capital_irfs.c_R_eps_nu(1)*(1-invest_share-cons_share_K)/cons_share

%%
% Basic TANK model, experiment with different parameters
sigma = 1.0;
phi=1.0;
phi_pi = 1.5;
phi_y  = 0.0;%.5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.97;
alpha=0.25;
epsilon=6;
lambda = 0.2;
Lambda = 0.0;
cons_share_K = (1-Lambda*(1-beta))*lambda*(1-alpha)*(epsilon-1)/epsilon;
cons_share_R = 1-cons_share_K;
labor_share_K = lambda;
labor_share_R = 1-labor_share_K;
% Loop over parameters
dynare TANKmodel_dgwages.mod noclearall;
options_.nograph = 1;
%param_loop = 0.0:0.05:2;

%% run different experiments
param_to_vary = [string('lambda'),string('Lambda')];
plot_xlabel = [string('Proportion of Keynesian Households'),string('Keynesian Household Debt Level')];
param_loop.(char(param_to_vary(1))) = 0.0:0.01:0.3;
param_loop.(char(param_to_vary(2))) = 0.0:0.05:2;
for j=1:length(param_to_vary)
    set_param_value('lambda',0.2);
    set_param_value('Lambda',0.0);
    set_param_value('sigma',2.0);
    param_array = param_loop.(char(param_to_vary(j)));
    Transmission_Channels = zeros(length(param_loop),7);
    suff_stats = zeros(length(param_array),4);
    YRP_changes = zeros(length(param_array),3);
    checks = zeros(length(param_array),3);
    for i=1:length(param_array)
        set_param_value(char(param_to_vary(j)),param_array(i));
        run 'CalcTransmissionChannels_TANKmodel.m';
        Transmission_Channels(i,1) = agg_inc*nominal_i_scale;
        Transmission_Channels(i,2) = het_inc*nominal_i_scale;
        Transmission_Channels(i,3) = ire*nominal_i_scale;
        Transmission_Channels(i,4) = fisher*nominal_i_scale;
        Transmission_Channels(i,5) = ies*nominal_i_scale;
        Transmission_Channels(i,6) = dC_C_TANK*nominal_i_scale;
        Transmission_Channels(i,7) = param_array(i);
        suff_stats(i,1) = Inc_wt_MPC_TANK;
        suff_stats(i,2) = Elas_R_TANK;
        suff_stats(i,3) = Elas_P_TANK;
        suff_stats(i,4) = Elas_EIS_TANK;
        YRP_changes(i,1) = dY_Y_TANK*nominal_i_scale;
        YRP_changes(i,2) = dR_R_TANK*nominal_i_scale;
        YRP_changes(i,3) = dP_P_TANK*nominal_i_scale;
        checks(i,1) = error_TANK;
        checks(i,2) = check_R;
        checks(i,3) = check_K;
    end
    cum_trans_channels = cumsum(Transmission_Channels(:,1:5),2);
    figure;
    area(param_array,Transmission_Channels(:,1:5));
    colormap(linspecer);
    legend('Agg Income', 'Het Income','Int. Rate Exposure', 'Fisher', 'Intertemporal Subs.','location','NorthWest');
    legend('boxoff');
    xlabel(char(plot_xlabel(j)));
    ylabel('Consumption Change');
    title('Transmission Channels to a 1% Interest Rate Decline');
end

