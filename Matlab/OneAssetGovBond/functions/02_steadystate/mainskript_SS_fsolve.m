
%% Set grids
[grid]  = makegrid(mpar,grid);

%% Use Tauchen method to approximate state space
[P_H,grid,par]=stochastics_variance(par, mpar,grid);

%% Solve for steady state

[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);
[c_guess,m_star,joint_distr,W_fc,Profits_fc,Output,N,grid,excess,par,C_agg,C_ind,inc,X_agg]...
    =steadystate_fsolve(P_H,grid,mpar,par,meshes);


%% Calculate steady state capital and further statistics
SS_stats

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