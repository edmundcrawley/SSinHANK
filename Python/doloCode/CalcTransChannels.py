"""
Calculates Auclert statistics for the TANK model
"""
import numpy as np
from dolo import *  # load the dolo library

def CalcTransChannels_TANK(model, T=8):
    #get parameters needed for later calculations
    param_names = model.symbols['parameters']
    param_vals = model.calibration['parameters']
    beta = param_vals[param_names.index('beta')]
    sigma = param_vals[param_names.index('sigma')]
    cons_share_K = param_vals[param_names.index('cons_share_K')]
    cons_share_R = param_vals[param_names.index('cons_share_R')]
    debt_limit = param_vals[param_names.index('debt_limit')]   
    #Solve model by pertubation
    dr_pert = approximate_controls(model, order=1)
    irf = response(model,dr_pert, 'eps_nu', T)
    #Now calc Auclert's statistics
    MPC_TANK_K = 1.0
    MPC_TANK_R = 1-beta
    y_K_share = 1/(1-debt_limit*(1-beta))*cons_share_K #share of total income earned by Keynesian labor
    y_R_share = 1- y_K_share
    MPC_TANK = MPC_TANK_R*y_R_share + MPC_TANK_K*y_K_share #income weighted MPC
    # 1) Aggregate income elasticity
    Inc_wt_MPC_TANK = MPC_TANK

    # 3) Interest Rate Exposure Channel (URE measures as a fraction of total
    # consumption/income)
    URE_K = -beta*debt_limit*y_K_share  #Nominal debt of Keynesians is a multiple of their labor income
    URE_R = - URE_K
    Elas_R_TANK = URE_K*MPC_TANK_K + URE_R*MPC_TANK_R
    # 4) Fisher channel (debt deflation)
    NNP_K = -debt_limit*y_K_share
    NNP_R = - NNP_K
    Elas_P_TANK = NNP_K*MPC_TANK_K + NNP_R*MPC_TANK_R
    # 5) Intertemporal Substitution Channel
    Hicks_scaling_TANK = (1.0-MPC_TANK_R)*cons_share_R + (1.0-MPC_TANK_K)*cons_share_K
    Elas_EIS_TANK = Hicks_scaling_TANK/sigma

    # Get inputs for partial eq. decomposition
    # dY/Y
    dY_Y_TANK = float(irf.sel(V='y_gap')[2])
    dYK_Y_TANK = (float(irf.sel(V='w_real')[2]) + float(irf.sel(V='n_K')[2]))*y_K_share
    dYR_Y_TANK = dY_Y_TANK - dYK_Y_TANK
    # dR/R
    dR_R_TANK = float(irf.sel(V='r_real')[2])
    # dP/P
    dP_P_TANK = float(irf.sel(V='pi')[2])

    # Get output for partial eq. decomposition
    # dC/C
    dC_C_TANK = float(irf.sel(V='y_gap')[2])

    # Check they add up...
    dC_C_Auclert_TANK = MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK \
                    + Elas_R_TANK*dR_R_TANK - Elas_P_TANK*dP_P_TANK \
                    - Elas_EIS_TANK*dR_R_TANK
    error_TANK = dC_C_Auclert_TANK - dC_C_TANK

    # Agg income channel
    agg_inc = MPC_TANK*dY_Y_TANK
    # Heterogeneous Income Channel
    het_inc = MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK - MPC_TANK*dY_Y_TANK
    # Unhedged Interest Rate Exposure
    ire = Elas_R_TANK*dR_R_TANK;
    # Fisher Channel
    fisher = - Elas_P_TANK*dP_P_TANK;
    # Intertemporal Elasticity Channel
    ies = - Elas_EIS_TANK*dR_R_TANK;

    #check for Keynesians
    dC_K = MPC_TANK_K*dYK_Y_TANK + URE_K*MPC_TANK_K*dR_R_TANK \
        - NNP_K*MPC_TANK_K*dP_P_TANK
    check_K = dC_K- float(irf.sel(V='c_K')[2])*cons_share_K

    #check for Ricardians
    dC_R = MPC_TANK_R*dYR_Y_TANK + URE_R*MPC_TANK_R*dR_R_TANK \
        - NNP_R*MPC_TANK_R*dP_P_TANK \
        - Elas_EIS_TANK*dR_R_TANK
    check_R = dC_R- float(irf.sel(V='c_R')[2])*(1-cons_share_K)

    #Scale for a 1% decrease in interest rates
    nominal_i_scale = -1.0/float(irf.sel(V='i')[2])
    
    Transmission_Channels = [agg_inc*nominal_i_scale, het_inc*nominal_i_scale, \
                             ire*nominal_i_scale, fisher*nominal_i_scale, \
                             ies*nominal_i_scale, dC_C_TANK*nominal_i_scale]
    suff_stats = [Inc_wt_MPC_TANK, Elas_R_TANK, Elas_P_TANK, Elas_EIS_TANK]
    YRP_changes = [dY_Y_TANK*nominal_i_scale, dR_R_TANK*nominal_i_scale, dP_P_TANK*nominal_i_scale]
    checks = [error_TANK, check_R, check_K]
    
    return Transmission_Channels, suff_stats, YRP_changes, checks
    
 

def CalcTransChannels_TANKcapital(model, T=8):
    #get parameters needed for later calculations
    param_names = model.symbols['parameters']
    param_vals = model.calibration['parameters']
    beta = param_vals[param_names.index('beta')]
    sigma = param_vals[param_names.index('sigma')]
    cons_share_K = param_vals[param_names.index('cons_share_K')]
    cons_share_R = param_vals[param_names.index('cons_share_R')]
    debt_limit = param_vals[param_names.index('debt_limit')]   
    invest_share = param_vals[param_names.index('invest_share')]   
    #Solve model by pertubation
    dr_pert = approximate_controls(model, order=1)
    irf = response(model,dr_pert, 'eps_nu', T)
    #Now calc Auclert's statistics
    MPC_TANK_K = 1.0
    MPC_TANK_R = 1-beta
    y_K_share = 1/(1-debt_limit*(1-beta))*cons_share_K #share of total income earned by Keynesian labor
    y_R_share = 1- y_K_share
    MPC_TANK = MPC_TANK_R*y_R_share + MPC_TANK_K*y_K_share #income weighted MPC
    # 1) Aggregate income elasticity
    Inc_wt_MPC_TANK = MPC_TANK

    # 3) Interest Rate Exposure Channel (URE measures as a fraction of total
    # consumption/income)
    cons_share = cons_share_R+cons_share_K
    URE_K = -beta*debt_limit*y_K_share  #Nominal debt of Keynesians is a multiple of their labor income
    URE_R = - URE_K + invest_share
    Elas_R_TANK = (URE_K*MPC_TANK_K + URE_R*MPC_TANK_R)/cons_share
    # 4) Fisher channel (debt deflation)
    NNP_K = -debt_limit*y_K_share
    NNP_R = - NNP_K
    Elas_P_TANK = (NNP_K*MPC_TANK_K + NNP_R*MPC_TANK_R)/cons_share
    # 5) Intertemporal Substitution Channel
    Hicks_scaling_TANK = (1.0-MPC_TANK_R)*cons_share_R/cons_share + (1.0-MPC_TANK_K)*cons_share_K/cons_share
    Elas_EIS_TANK = Hicks_scaling_TANK/sigma

    # Get inputs for partial eq. decomposition
    # dY/Y
    dY_Y_TANK = float(irf.sel(V='y')[2])
    dYK_Y_TANK = (float(irf.sel(V='w_real')[2]) + float(irf.sel(V='n_K')[2]))*y_K_share
    dYR_Y_TANK = dY_Y_TANK - dYK_Y_TANK
    # dR/R
    dR_R_TANK = float(irf.sel(V='r_real')[2])
    # dP/P
    dP_P_TANK = float(irf.sel(V='pi')[2])

    # Get output for partial eq. decomposition
    # dC/C
    dC_C_TANK = float(irf.sel(V='c_R')[2])*cons_share_R/cons_share + float(irf.sel(V='c_K')[2])*cons_share_K/cons_share
    dC_C_TANK_R = float(irf.sel(V='c_R')[2])*cons_share_R/cons_share_R
    dC_C_TANK_K = float(irf.sel(V='c_K')[2])*cons_share_K/cons_share_K

    # Check they add up...
    dC_C_Auclert_TANK = (MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK)/cons_share \
                    + Elas_R_TANK*dR_R_TANK - Elas_P_TANK*dP_P_TANK \
                    - Elas_EIS_TANK*dR_R_TANK
    error_TANK = dC_C_Auclert_TANK - dC_C_TANK

    # Agg income channel
    agg_inc = MPC_TANK*dY_Y_TANK/cons_share
    # Heterogeneous Income Channel
    het_inc = (MPC_TANK_R*dYR_Y_TANK + MPC_TANK_K*dYK_Y_TANK - MPC_TANK*dY_Y_TANK)/cons_share
    # Unhedged Interest Rate Exposure
    ire = Elas_R_TANK*dR_R_TANK;
    # Fisher Channel
    fisher = - Elas_P_TANK*dP_P_TANK;
    # Intertemporal Elasticity Channel
    ies = - Elas_EIS_TANK*dR_R_TANK;

    #check for Keynesians
    dC_K = MPC_TANK_K*dYK_Y_TANK/cons_share \
        + URE_K*MPC_TANK_K*dR_R_TANK/cons_share \
        - NNP_K*MPC_TANK_K*dP_P_TANK/cons_share
    check_K = dC_K- float(irf.sel(V='c_K')[2])*cons_share_K/cons_share

    #check for Ricardians
    dC_R = MPC_TANK_R*dYR_Y_TANK/cons_share \
        + URE_R*MPC_TANK_R*dR_R_TANK/cons_share \
        - NNP_R*MPC_TANK_R*dP_P_TANK/cons_share \
        - Elas_EIS_TANK*dR_R_TANK
    check_R = dC_R- float(irf.sel(V='c_R')[2])*(1-invest_share-cons_share_K)/cons_share

    #Scale for a 1% decrease in interest rates
    nominal_i_scale = -1.0/float(irf.sel(V='i')[2])
    
    Transmission_Channels = [agg_inc*nominal_i_scale, het_inc*nominal_i_scale, \
                             ire*nominal_i_scale, fisher*nominal_i_scale, \
                             ies*nominal_i_scale, dC_C_TANK*nominal_i_scale]
    suff_stats = [Inc_wt_MPC_TANK, Elas_R_TANK, Elas_P_TANK, Elas_EIS_TANK]
    YRP_changes = [dY_Y_TANK*nominal_i_scale, dR_R_TANK*nominal_i_scale, dP_P_TANK*nominal_i_scale]
    checks = 100*np.array([error_TANK/dC_C_TANK, check_R/dC_C_TANK_R, check_K/dC_C_TANK_K])
    
    IRF_i = np.array(irf.sel(V='i')[2:])*nominal_i_scale
    IRF_c_R = np.array(irf.sel(V='c_R')[2:])*nominal_i_scale
    IRF_c_K = np.array(irf.sel(V='c_K')[2:])*nominal_i_scale
    IRF_k = np.array(irf.sel(V='k')[2:])*nominal_i_scale
    IRF_q = np.array(irf.sel(V='q')[2:])*nominal_i_scale
    IRF_r_real = np.array(irf.sel(V='r_real')[2:])*nominal_i_scale
    
    return Transmission_Channels, suff_stats, YRP_changes, checks, IRF_i, IRF_c_R, IRF_c_K, IRF_k, IRF_q, IRF_r_real
    