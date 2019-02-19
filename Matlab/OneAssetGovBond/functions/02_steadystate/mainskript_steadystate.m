
%% Set grids
[grid]  = makegrid(mpar,grid);

%% Use Tauchen method to approximate state space
[P_H,grid,par]=stochastics_variance(par, mpar,grid);

%% Solve for steady state
rmax = 1/par.beta-1 ; % max is complete market r
rmin = -0.03;
r = (rmax+rmin)/2;      % initial r for bisection
par.RB=1+r;
init = 999;             % initial gap between two r for bisection

[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

while abs(init) > mpar.crit
    
    [ N,W_fc,Profits_fc,WW,RBRB,Output] = factor_returns(meshes,grid,par,mpar);
    
    [c_guess,inc]=policyguess(meshes,WW,RBRB,par,mpar);
    count=0;
    
    % Solve Policies and Joint Distribution
    disp('Solving household problem by EGM')
    tic
    [c_guess,m_star,distPOL]=...
        policies_SS(c_guess, grid, inc,RBRB,P_H,mpar,par);
    toc
    
    disp(([distPOL]));
    
    disp('Calc Joint Distr')
    
    [joint_distr]=JDiteration(m_star,P_H,mpar,grid);
    joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
    
    aux_x = par.tau*N/par.H*W_fc*meshes.h/(1+par.gamma);
    aux_x(:,end)=0;
    C_ind=c_guess+aux_x;
    C_agg = joint_distr(:)'*C_ind(:);
    
    AggregateSavings=m_star(:)'*joint_distr(:);
    ExcessA=AggregateSavings;
    
    % Use Bisection to update r
    if ExcessA > 0
        rmax = (r+rmax)/2;
    else
        rmin = (r+rmin)/2;
    end
    init = rmax-rmin;
    disp('Starting Iteration for r. Difference remaining:                      ');
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %20.9f \n',init);
    r = (rmax+rmin)/2;
    par.RB = 1+r;
    
end

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