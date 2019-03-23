%%%%%%%%%%%%%% Vary over the parameter psi
set_param_value('psi_c',psi_c);
%%%%%%%%%%%%%%%

stoch_simul(var_list_);
TANK_Capital_irfs = oo_.irfs;

% Now calc Auclert's statistics
MPC_TANK_Capital_K = 1.0;
MPC_TANK_Capital_R = 1-beta;
y_K_share = 1/(1-Lambda*(1-beta))*cons_share_K; %share of total income earned by Keynesian labor
y_R_share = 1- y_K_share;
MPC_TANK_Capital = MPC_TANK_Capital_R*y_R_share + MPC_TANK_Capital_K*y_K_share;
% 1) Aggregate income elasticity
Inc_wt_MPC_TANK_Capital = MPC_TANK_Capital;
% 2) Income Heterogeneity Channel
% ????
% 3) Interest Rate Exposure Channel (URE measures as a fraction of total
% consumption/income)
cons_share = cons_share_R+cons_share_K;
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
dC_C_TANK_Capital_R = TANK_Capital_irfs.c_R_eps_nu(1)*cons_share_R/cons_share_R;
dC_C_TANK_Capital_K = TANK_Capital_irfs.c_K_eps_nu(1)*cons_share_K/cons_share_K;


% Check they add up...
dC_C_Auclert_TANK_Capital = (MPC_TANK_Capital_R*dYR_Y_TANK_Capital + MPC_TANK_Capital_K*dYK_Y_TANK_Capital)/cons_share ...
                    + Elas_R_TANK_Capital*dR_R_TANK_Capital - Elas_P_TANK_Capital*dP_P_TANK_Capital ... 
                    - Elas_EIS_TANK_Capital*dR_R_TANK_Capital;
error_TANK_Capital = dC_C_Auclert_TANK_Capital - dC_C_TANK_Capital;

% Agg income channel
agg_inc = (MPC_TANK_Capital*dY_Y_TANK_Capital)/cons_share;
% Heterogeneous Income Channel
het_inc = (MPC_TANK_Capital_R*dYR_Y_TANK_Capital + MPC_TANK_Capital_K*dYK_Y_TANK_Capital - MPC_TANK_Capital*dY_Y_TANK_Capital)/cons_share;
% Unhedged Interest Rate Exposure
ire = Elas_R_TANK_Capital*dR_R_TANK_Capital;
% Fisher Channel
fisher = - Elas_P_TANK_Capital*dP_P_TANK_Capital;
% Intertemporal Elasticity Channel
ies = - Elas_EIS_TANK_Capital*dR_R_TANK_Capital;

%check for Keynesians
dC_K = MPC_TANK_Capital_K*dYK_Y_TANK_Capital/cons_share ...
        + URE_K*MPC_TANK_Capital_K*dR_R_TANK_Capital/cons_share ...
        - NNP_K*MPC_TANK_Capital_K*dP_P_TANK_Capital/cons_share;
check_K = dC_K- TANK_Capital_irfs.c_K_eps_nu(1)*cons_share_K/cons_share;

%check for Ricardians
dC_R = MPC_TANK_Capital_R*dYR_Y_TANK_Capital/cons_share ...
        + URE_R*MPC_TANK_Capital_R*dR_R_TANK_Capital/cons_share ...
        - NNP_R*MPC_TANK_Capital_R*dP_P_TANK_Capital/cons_share ...
        - Elas_EIS_TANK_Capital*dR_R_TANK_Capital;
check_R = dC_R- TANK_Capital_irfs.c_R_eps_nu(1)*(1-invest_share-cons_share_K)/cons_share;

%Scale for a 1% decrease in interest rates
nominal_i_scale = -1/TANK_Capital_irfs.i_eps_nu(1);
