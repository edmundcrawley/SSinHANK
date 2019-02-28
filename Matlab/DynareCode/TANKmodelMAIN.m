%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for TANK model in Dynare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('c:\dynare'))

dynare 'RepAgent.mod'

% MPC in this representative agent model is (1-beta)/(1+sigma/phi*(1-alpha))
MPC_RANK = (1-beta)/(1+sigma/phi*(1-alpha));
% Size of two channels (no other channels in rep agent model)
% 1) Aggregate income elasticity
Inc_wt_MPC_RANK = MPC_RANK;
% 2) Intertemporal Substitution Channel
Hicks_scaling_RANK = 1.0-MPC_RANK;
Elas_R_RANK = Hicks_scaling_RANK*sigma;

% Get inputs for partial eq. decomposition
% dY/Y
dY_Y = oo_.irfs.y_gap_eps_nu(1);
% dR/R
dR_R = oo_.irfs.r_real_eps_nu(1);

% Get output for partial eq. decomposition
% dC/C
dC_C = oo_.irfs.y_gap_eps_nu(1);

% Check they add up...
dC_C_Auclert = Inc_wt_MPC_RANK*dY_Y - Elas_R_RANK*dR_R;
error = dC_C_Auclert - dC_C;


