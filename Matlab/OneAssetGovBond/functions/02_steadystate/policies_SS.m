function [c_new,m_star,distPOL] = policies_SS(c_guess, grid, inc,RBRB,P,mpar,par)

%% Apply EGM to solve for optimal policies and marginal utilities
money_expense  = repmat(grid.m',[1 mpar.nh]);
distC = 99999;

count=0;
while max([distC])>mpar.crit
    count=count+1;
    
    %% Update policies 
    mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy 
    
    aux=reshape(permute(mutil_c,[2 1]),[mpar.nh mpar.nm]);
    % form expectations
    EMU_aux = par.beta*RBRB.*permute(reshape(P*aux,[mpar.nh mpar.nm]),[2 1]);
    
    c_aux = 1./(EMU_aux.^(1/par.xi));
    
    % Take budget constraint into account
    [c_new,m_star]=EGM(grid,inc,money_expense,c_aux,mpar,par);
    
    m_star(m_star>grid.m(end)) = grid.m(end); %no extrapolation
    
    %% Step 6: Check convergence of policies
    distC = max((abs(c_guess(:)-c_new(:))));
    
    % Update c policy guesses
    c_guess=c_new;
    
end
distPOL=[distC];
end

%% SUBFUNCTIONS

function [c_update,m_update]=EGM(grid,inc,money_expense,c_aux,mpar,par)
%%EGM_Step1_b computes the optimal consumption and corresponding optimal bond holdings by taking the budget constraint into account.
% c_update(m,h):    Update for consumption policy 
% m_update(m,h):    Update for bond policy 

%% EGM: Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_star = (c_aux + money_expense - inc.labor - inc.profits);
RR = (par.RB+(m_star<0)*par.borrwedge)./par.PI;
m_star = m_star./RR;

% Identify binding constraints
binding_constraints = money_expense < repmat(m_star(1,:),[mpar.nm 1 ]);

% Consumption when drawing assets m' to zero: Eat all Resources
Resource = inc.labor  + inc.money + inc.profits;

%% Next step: Interpolate on grid

c_update=zeros(mpar.nm,mpar.nh);
m_update=zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
    Savings=griddedInterpolant(m_star(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_update(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_star(:,hh),c_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_update(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end

c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_update(binding_constraints) = min(grid.m);

end
