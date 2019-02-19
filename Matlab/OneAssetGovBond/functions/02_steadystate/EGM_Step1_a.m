
function [c_update,m_update]=EGM_Step1_a(grid,inc,money_expense,c_aux,mpar,par)
%%EGM_Step1_b computes the optimal consumption and corresponding optimal money
% holdings in case the capital stock cannot be adjusted by taking the budget
% constraint into account.
% c_update(m,k,s*h,M,K):    Update for consumption policy under no-adj.
% m_update(m,k,s*h,M,K):    Update for money policy under no-adj.

%% EGM: Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star = (c_aux + money_expense - inc.labor - inc.profits);
RR = (par.RB+(m_star<0)*par.borrwedge)./par.PI;
m_star = m_star./RR;
% m_star is old saving

% Identify binding constraints
binding_constraints = money_expense < repmat(m_star(1,:),[mpar.nm 1 ]);
%                    -------------- 
%                     new saving

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor  + inc.money + inc.profits;

%% Next step: Interpolate w_guess and c_guess from new k-grids
% using c(s,h,k',K), k(s,h,k'Kf

m_star = reshape(m_star,[mpar.nm mpar.nh]);
c_aux= reshape(c_aux,[mpar.nm mpar.nh]);

%Interpolate grid.m and c_n_aux defined on m_star_n over grid.m
% [c_update, m_update]=egm1b_aux_mex(grid.m,m_star,c_aux);

c_update = zeros(mpar.nm,mpar.nh);
m_update = zeros(mpar.nm,mpar.nh);

for s=1:mpar.nh % produce consumption, saving policy function for each income grid (grid.s)
    M_choice=griddedInterpolant(m_star(:,s),grid.m); % generate M choice function m(m*,k,h)=m' 
    m_update(:,s)=M_choice(grid.m); % Obtain m'(h,k,m) by Interpolation (notice this is out of grid, used linear interpolation)
    % This m_update is the interpolated (true) policy function : M(h,m)
    C_choice=griddedInterpolant(m_star(:,s),c_aux(:,s)); % generate consumption function c(s,a*(s,a'))
    c_update(:,s)=C_choice(grid.m);  % Obtain c(s,a) by interpolation (notice this is out of grid, used linear interpolation)
    % This C is the interpolated (true) policy function : C(s,a)           
end

c_update = reshape(c_update,[mpar.nm, mpar.nh]);
m_update = reshape(m_update,[mpar.nm, mpar.nh]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new handling of no borrowing
binding_constraints = c_update>Resource;
c_update(binding_constraints) = Resource(binding_constraints);
m_update(binding_constraints) = min(grid.m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
%m_update(binding_constraints) = min(grid.m);

end