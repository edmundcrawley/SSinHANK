function [Difference,LHS,RHS,JD_new,c_star,m_star,P]  = Fsys_Gov_new(State,Stateminus,...
    Control_sparse,Controlminus_sparse,StateSS,...
    ControlSS,Gamma_state,Gamma_control,InvGamma,Copula,...
    par,mpar,grid,targets,P,aggrshock,oc)
% System of equations written in Schmitt-Grohé-Uribe generic form with states and controls
% STATE: Vector of state variables t+1 (only marginal distributions for histogram)
% STATEMINUS: Vector of state variables t (only marginal distributions for histogram)
% CONTROL: Vector of state variables t+1 (only coefficients of sparse polynomial)
% CONTROLMINUS: Vector of state variables t (only coefficients of sparse polynomial)
% STATESS and CONTROLSS: Value of the state and control variables in steady
% state. For the Value functions these are at full grids.
% GAMMA_STATE: Mapping such that perturbationof marginals are still
% distributions (sum to 1).
% PAR, MPAR: Model and numerical parameters (structure)
% GRID: Liquid, illiquid and productivity grid
% TARGETS: Stores targets for government policy
% P: steady state transition matrix
% aggrshock: sets wether the Aggregate shock is TFP or uncertainty

%% Initializations
mutil = @(c)(1./(c.^par.xi));
invmutil = @(mu)((1./mu).^(1/par.xi));
% Generate meshes for b,k,h
[meshes.m,meshes.h] = ndgrid(grid.m,grid.h);

% Number of states, controls
nx   = mpar.numstates; % Number of states
ny   = mpar.numcontrols; % number of Controls
NxNx = nx-2; % Number of states without aggregate
NN   = mpar.nm*mpar.nh; % Number of points in the full grid

% Initialize LHS and RHS
LHS  = zeros(nx+ny,1);
RHS  = zeros(nx+ny,1);

%% Indexes for LHS/RHS
% Indexes for Controls
mutil_cind = 1:NN;
PIind = 1*NN+1;
Yind  = 1*NN+2;
Wind  = 1*NN+3;
Profitind = 1*NN+4;
Nind  = 1*NN+5;
Bind  = 1*NN+6;
Cind  = 1*NN+7;
Gind  = 1*NN+8;
Xind  = 1*NN+9;
MCind  = 1*NN+10;

% Indexes for States
marginal_mind = (1:mpar.nm-1);
marginal_hind = (mpar.nm-1 + (1:(mpar.nh-2)));

RBind = NxNx+1;
Sind  = NxNx+2;

%% Control Variables (Change Value functions according to sparse polynomial)
Control      = ControlSS .* (1+Gamma_control*(Control_sparse));
Controlminus = ControlSS .* (1+Gamma_control*(Controlminus_sparse));

Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Control_sparse);
Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Controlminus_sparse);

%% State Variables
% read out marginal histogramm in t+1, t
Distribution      = StateSS(1:end-2) + Gamma_state * State(1:NxNx);
Distributionminus = StateSS(1:end-2) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Endogenous States
RB      = StateSS(end-1) + (State(end-1));
RBminus = StateSS(end-1) + (Stateminus(end-1));

% Aggregate Exogenous States
S       = StateSS(end) + (State(end));
Sminus  = StateSS(end) + (Stateminus(end));

%% Split the Control vector into items with names
% Controls
mutil_c       = mutil(Control(mutil_cind));
mutil_cminus  = mutil(Controlminus(mutil_cind));

% Aggregate Controls (t+1)
PI = exp(Control(PIind));
Y  = exp(Control(Yind ));
B  = exp(Control(Bind ));

% Aggregate Controls (t)
PIminus = exp(Controlminus(PIind));
Yminus  = exp(Controlminus(Yind ));
Wminus  = exp(Controlminus(Wind ));
Profitminus  = exp(Controlminus(Profitind ));
Nminus  = exp(Controlminus(Nind ));
Bminus  = exp(Controlminus(Bind ));
Cminus  = exp(Controlminus(Cind ));
Gminus  = Controlminus(Gind );
% Tminus  = exp(Controlminus(Tind ));
Xminus  = exp(Controlminus(Xind ));
MCminus  = exp(Controlminus(MCind ));
%% Write LHS values
% Controls
LHS(nx+mutil_cind) = invmutil(mutil_cminus);
LHS(nx+Yind)       = (Yminus);
LHS(nx+Wind)       = (Wminus);
LHS(nx+Profitind)  = (Profitminus);
LHS(nx+Nind)       = (Nminus);
LHS(nx+Bind)       = (Bminus);
LHS(nx+Cind)       = (Cminus);
LHS(nx+Gind)       = (Gminus);
% LHS(nx+Tind)       = (Tminus);
LHS(nx+Xind)       = (Xminus);
LHS(nx+MCind)       = (MCminus);
% States
% Marginal Distributions (Marginal histograms)
LHS(marginal_mind) = Distribution(1:mpar.nm-1);
LHS(marginal_hind) = Distribution(mpar.nm+(1:mpar.nh-2));

LHS(RBind)         = (RB);
LHS(Sind)          = (S);

% take into account that RB is in logs
RB=exp(RB); RBminus=exp(RBminus);

%% Set of Differences for exogenous process
RHS(Sind) = (par.rhoS * (Sminus));

switch(aggrshock)
    case('MP')
        EPS_TAYLOR=Sminus;
        TFP=1;
    case('TFP')
        TFP=exp(Sminus);
        EPS_TAYLOR=0;
    case('Uncertainty')
        TFP=1;
        EPS_TAYLOR=0;
        % Tauchen style for Probability distribution next period
        [P,~,~] = ExTransitions(exp(Sminus),grid,mpar,par);
end

marginal_mminus = Distributionminus(1:mpar.nm)';
marginal_hminus = Distributionminus(mpar.nm+(1:mpar.nh))';

Hminus  = sum(grid.h(1:end-1).*marginal_hminus(1:end-1)); %Last column is entrepreneurs.
Lminus  = sum(grid.m.*marginal_mminus);

RHS(nx+Bind) = Lminus;

% Calculate joint distributions
cumdist = zeros(mpar.nm+1,mpar.nh+1);
cumdist(2:end,2:end) = Copula({cumsum(marginal_mminus),cumsum(marginal_hminus)});
JDminus = diff(diff(cumdist,1,1),1,2);

%% Aggregate Output
RHS(nx+MCind)             =  par.mu- (par.beta * log(PI)*Y/Yminus - log(PIminus))/par.kappa;

RHS(nx+Nind)    =  (par.tau*TFP*par.alpha*grid.K.^(1-par.alpha).*MCminus).^(1/(1-par.alpha+par.gamma));

RHS(nx+Yind)    = (TFP*(Nminus).^(par.alpha).*grid.K.^(1-par.alpha));

% Wage Rate
RHS(nx+Wind)    = TFP *par.alpha*MCminus.* (grid.K./(Nminus)).^(1-par.alpha);

% Profits for Entrepreneurs
RHS(nx+Profitind) = (1-MCminus)*Yminus - Yminus.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2;

%% Wages net of leisure services
WW=par.gamma/(1+par.gamma)*(Nminus./Hminus)*Wminus*ones(mpar.nm,mpar.nh);
WW(:,end)=Profitminus*par.profitshare;

%% Incomes (grids)
inc.labor   = par.tau*WW.*(meshes.h);
inc.money   = meshes.m.*(RBminus/PIminus+(meshes.m<0).*par.borrwedge/PIminus);

jd_aux = sum(JDminus,1);
inc.profits = sum((1-par.tau)*par.gamma/(1+par.gamma).*(Nminus/par.H).*Wminus.*grid.h(1:end-1).*jd_aux(1:end-1)) ...
                + (1-par.tau)*Profitminus*par.profitshare*jd_aux(end) - (RBminus/PIminus*Bminus-B); % lump sum transfer

%% Update policies
RBaux=(RB+(meshes.m<0).*par.borrwedge)/PI;
EVm = reshape(reshape(RBaux(:).*mutil_c,[mpar.nm mpar.nh])*P',[mpar.nm, mpar.nh]);

[c_star,m_star] = EGM_policyupdate(EVm,PIminus,RBminus,inc,meshes,grid,par,mpar);

aux_x=par.tau*(Nminus/Hminus)*Wminus*meshes.h / (1+par.gamma);
% aux_x: % consumption in return for labor supply
aux_x(:,end)=0; % Ent's consumption in return for labor supply = 0

c=c_star+aux_x; % c_star: net consumption (on composite goods)
RHS(nx+Cind) = JDminus(:)'*c(:);
RHS(nx+Xind) = JDminus(:)'*c_star(:);

%% Update Marginal Value Bonds
mutil_c_aux = mutil(c_star); % marginal utility at consumption policy no adjustment
RHS(nx+mutil_cind) = invmutil(mutil_c_aux(:)); % Write Marginal Utility to RHS of F

%% Differences for distributions
% find next smallest on-grid value for money choices
weight11  = zeros(mpar.nm, mpar.nh,mpar.nh);
weight12  = zeros(mpar.nm, mpar.nh,mpar.nh);

% Adjustment case
[Dist_m,idm] = genweight(m_star,grid.m);

idm=repmat(idm(:),[1 mpar.nh]);
idh=kron(1:mpar.nh,ones(1,mpar.nm*mpar.nh));

index11 = sub2ind([mpar.nm mpar.nh],idm(:),idh(:));
index12 = sub2ind([mpar.nm mpar.nh],idm(:)+1,idh(:));

for hh=1:mpar.nh
    
    %Corresponding weights
    weight11_aux = (1-Dist_m(:,hh));
    weight12_aux =  (Dist_m(:,hh));
    
    % Dimensions (mxk,h',h)
    weight11(:,:,hh)=weight11_aux(:)*P(hh,:);
    weight12(:,:,hh)=weight12_aux(:)*P(hh,:);
end

weight11=permute(weight11,[1 3 2]);
weight12=permute(weight12,[1 3 2]);

rowindex=repmat(1:mpar.nm*mpar.nh,[1 2*mpar.nh]);

H=sparse(rowindex,[index11(:); index12(:)],...
    [weight11(:); weight12(:)],mpar.nm*mpar.nh,mpar.nm*mpar.nh); % mu'(h',k'), a without interest


JD_new=JDminus(:)'*H;
JD_new = reshape(JD_new(:),[mpar.nm,mpar.nh]);
% Next period marginal histograms
% liquid assets
aux_m = squeeze(sum(JD_new,2));
RHS(marginal_mind) = aux_m(1:end-1); %Leave out last state
% human capital
aux_h = squeeze(sum(JD_new,1));
RHS(marginal_hind) = aux_h(1:end-2); %Leave out last state & entrepreneurs

%% Third Set: Government Budget constraint
% Return on bonds (Taylor Rule)
RHS(RBind) = log(par.RB) + par.rho_R* log(RBminus/par.RB)  + log(PIminus/par.PI).*((1-par.rho_R)*par.theta_pi) ...
             + EPS_TAYLOR;

% % Inflation jumps to equilibrate real bond supply and demand
LHS(nx+PIind) = par.RB*targets.B;
RHS(nx+PIind) = RB/PI * B;


% LHS(nx+PIind) = log((B)/(targets.B));
% 
% RHS(nx+PIind) = par.rho_B * log((Bminus)/(targets.B)) ...
%     + par.rho_B * log(RBminus/par.RB)...
%     - (par.rho_B+par.gamma_pi) * log(PIminus/par.PI); % ... - par.gamma_T * log((Tminus)/(targets.T));

% Government expenditures
% RHS(nx+Gind) =  B - Bminus*RBminus/PIminus + Tminus;
RHS(nx+Gind) =  B - Bminus*RBminus/PIminus;

% taxrevenue =(1-par.tau).*Wminus.*Nminus +(1-par.tau).*Profitminus;

% RHS(nx+Tind) = taxrevenue;
%% Difference
Difference = InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);(ControlSS(1:end-oc));ones(oc,1)]);

end
function [ weight,index ] = genweight( x,xgrid )
% function: GENWEIGHT generates weights and indexes used for linear interpolation
%
% X: Points at which function is to be interpolated.
% xgrid: grid points at which function is measured
% no extrapolation allowed
[~,index] = histc(x,xgrid);
index(x<=xgrid(1))=1;
index(x>=xgrid(end))=length(xgrid)-1;

weight = (x-xgrid(index))./(xgrid(index+1)-xgrid(index)); % weight xm of higher gridpoint
weight(weight<=0) = 1.e-16; % no extrapolation
weight(weight>=1) = 1-1.e-16; % no extrapolation

end  % function
