

% % Calc steady state share of labor and consumption of each type
% cons_share_to_labor_share_K = (1-Lambda*(1-beta))*(epsilon-1)/epsilon*(1-alpha);
% cons_share_K_obj = @(x)x^sigma * (x/cons_share_to_labor_share_K)^phi - (lambda/(1-lambda))^(sigma+phi)*(1-x)^sigma * (1-x/cons_share_to_labor_share_K)^phi;
% cons_share_K = fsolve(cons_share_K_obj, lambda);
% %cons_share_K =cons_share_K(1);
% cons_share_R = 1-cons_share_K;
% labor_share_K = cons_share_K/cons_share_to_labor_share_K;
% labor_share_R = 1-labor_share_K;
% 
% set_param_value('cons_share_K',cons_share_K);
% set_param_value('cons_share_R',cons_share_R);
% set_param_value('labor_share_K',labor_share_K);
% set_param_value('labor_share_R',labor_share_R);
%%%%%%%%%%%%%%
cons_share_K = (1-Lambda*(1-beta))*lambda*(1-alpha)*(epsilon-1)/epsilon;
cons_share_R = 1-cons_share_K;
labor_share_K = lambda;
labor_share_R = 1-labor_share_K;

set_param_value('cons_share_K',cons_share_K);
set_param_value('cons_share_R',cons_share_R);
set_param_value('labor_share_K',labor_share_K);
set_param_value('labor_share_R',labor_share_R);
% %%%%%%%%%%%%%%


stoch_simul(var_list_);
%dynare 'TANKmodel.mod' noclearall;
TANK_irfs = oo_.irfs;

% Now calc Auclert's statistics
MPC_TANK_K = 1.0;
MPC_TANK_R = 1-beta;
y_K_share = 1/(1-Lambda*(1-beta))*cons_share_K; %share of total income earned by Keynesian labor
y_R_share = 1- y_K_share;
MPC_TANK = MPC_TANK_R*y_R_share + MPC_TANK_K*y_K_share; %income weighted MPC
% 1) Aggregate income elasticity
Inc_wt_MPC_TANK = MPC_TANK;
% 2) Income Heterogeneity Channel
% Calculate this exactly (rather than finding gamma)
y_K = TANK_irfs.w_real_eps_nu(1) + TANK_irfs.n_K_eps_nu(1);
y_R = TANK_irfs.w_real_eps_nu(1) + TANK_irfs.n_R_eps_nu(1);
IncomeChannel_Total = y_K*cons_share_K + y_R*(1-beta)*cons_share_R;
% 3) Interest Rate Exposure Channel (URE measures as a fraction of total
% consumption/income)
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
agg_inc = MPC_TANK*dY_Y_TANK;
% Heterogeneous Income Channel
het_inc = MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK - MPC_TANK*dY_Y_TANK;
% Unhedged Interest Rate Exposure
ire = Elas_R_TANK*dR_R_TANK;
% Fisher Channel
fisher = - Elas_P_TANK*dP_P_TANK;
% Intertemporal Elasticity Channel
ies = - Elas_EIS_TANK*dR_R_TANK;

%check for Keynesians
dC_K = MPC_TANK_K*dYK_Y_TANK + URE_K*MPC_TANK_K*dR_R_TANK ...
        - NNP_K*MPC_TANK_K*dP_P_TANK;
check_K = dC_K- TANK_irfs.c_K_eps_nu(1)*cons_share_K;

%check for Ricardians
dC_R = MPC_TANK_R*dYR_Y_TANK + URE_R*MPC_TANK_R*dR_R_TANK ...
        - NNP_R*MPC_TANK_R*dP_P_TANK ...
        - Elas_EIS_TANK*dR_R_TANK;
check_R = dC_R- TANK_irfs.c_R_eps_nu(1)*(1-cons_share_K);

%Scale for a 1% decrease in interest rates
nominal_i_scale = -1/TANK_irfs.i_eps_nu(1);