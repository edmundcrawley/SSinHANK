%% Parameters
% Household Parameters
par.beta        = 0.985;     % Discount factor
par.xi          = 1;          % CRRA
par.gamma       = 1000;          % Inverse Frisch elasticity

% Income Process
par.rhoH        = 0.979;    % Persistence of productivity
par.sigmaH      = 0.059;    % STD of productivity shocks
mpar.in         = 0.00046;  % Prob. to become entrepreneur
mpar.out        = 0.0625;   % Prob. to become worker again

% Firm Side Parameters
par.eta         = 20;
par.mu          = (par.eta-1)/par.eta;       % Markup 5%
par.alpha       = 1;%2/3/par.mu;  % Labor share 2/3

% Phillips Curve
par.prob_priceadj = 3/4; % average price duration of 4 quarters = 1/(1-par.prob_priceadj)
par.kappa         = (1-par.prob_priceadj)*(1-par.prob_priceadj*par.beta)/par.prob_priceadj;% Phillips-curve parameter (from Calvo prob.)

% Central Bank Policy
par.theta_pi    = 1.25; % Reaction to inflation
par.rho_R       = 0.5;  % Inertia


% Tax Schedule
par.tau         = 0.99;   % Proportional tax on labor and profit income 

% Debt rule
par.gamma_pi    = 1.5;   % 1.5, Reaction to inflation
par.gamma_T     = 0.5075; % Reaction to tax revenue
par.rho_B       = 0.99;  % Autocorrelation

%% Returns
par.PI          = 1;          % Gross inflation
par.RB          = 1;          % Market clearing interest rate to be solved for
par.borrwedge   = 1.02^0.25-1; % Wedge on borrowing

%% Grids
% Idiosyncratic States
mpar.nm         = 100;
mpar.nh         = 4;
mpar.tauchen    ='importance';

grid.K = 1; 

%% Numerical Parameters
mpar.crit    = 1e-11;
