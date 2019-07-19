function [Difference,LHS,RHS]  = Fsys_RANK_SGU_order(State,Stateminus,Control_sparse,Controlminus_sparse,...
    StateSS,ControlSS,Gamma_control,InvGamma,par,mpar,grid,aggrshock,oc,os)
% ControlSS=Yss;StateSS=Xss;
% Stateminus=State_m; 
% Control_sparse=Contr;Controlminus_sparse=Contr_m;
%% Initializations

% Number of states, controls
nx   = mpar.numstates; % Number of states
ny   = mpar.numcontrols; % number of Controls

% Initialize LHS and RHS
LHS  = zeros(nx+ny,1);
RHS  = zeros(nx+ny,1);

%% Indexes for LHS/RHS
% Indexes for States
NxNx=0;
Kind = NxNx+1;
Qind = NxNx+2;
VRind = NxNx+3;
NWind = NxNx+4;
PHIMind = NxNx+5;
RBind = NxNx+6;
RRBind = NxNx+7;
TFPind = NxNx+8;
EMind = NxNx+9;
ENind = NxNx+10; 

% Indexes for Controls
NN=0;
YMind = 1*NN+1;
Nind  = 1*NN+2;
Iind = 1*NN+3;
Cind = 1*NN+4;
NWEind = 1*NN+5;
NWNind = 1*NN+6;
MCind = 1*NN+7;
Wind  = 1*NN+8;
Profitind  = 1*NN+9;
Yind = 1*NN+10;
LDind = 1*NN+11;
RBBind = 1*NN+12;
NUind = 1*NN+13;
ETAind = 1*NN+14;
Zind = 1*NN+15;
Xind = 1*NN+16;
PIind = 1*NN+17;

%% Control Variables (Change Value functions according to sparse polynomial)
% Control      = ControlSS .* (1+Gamma_control*(Control_sparse));
% Controlminus = ControlSS .* (1+Gamma_control*(Controlminus_sparse));
% Control(end-oc+1:end)       = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Control_sparse);
% Controlminus(end-oc+1:end)  = ControlSS(end-oc+1:end) + Gamma_control(end-oc+1:end,:)*(Controlminus_sparse);

% Control       = ControlSS + Gamma_control*(Control_sparse);
% Controlminus  = ControlSS + Gamma_control*(Controlminus_sparse);

Control       = ControlSS + (Control_sparse);
Controlminus  = ControlSS + (Controlminus_sparse);

%% State Variables
% read out marginal histogramm in t+1, t
% Distribution      = StateSS(1:end-os) + Gamma_state * State(1:NxNx);
% Distributionminus = StateSS(1:end-os) + Gamma_state * Stateminus(1:NxNx);

% Aggregate Endogenous States
K      = StateSS(end-9) + (State(end-9));
Kminus = StateSS(end-9) + (Stateminus(end-9));

Q      = StateSS(end-8) + (State(end-8));
Qminus = StateSS(end-8) + (Stateminus(end-8));

VR      = StateSS(end-7) + (State(end-7));
VRminus = StateSS(end-7) + (Stateminus(end-7));

NW      = StateSS(end-6) + (State(end-6));
NWminus = StateSS(end-6) + (Stateminus(end-6));

PHIM      = StateSS(end-5) + (State(end-5));
PHIMminus = StateSS(end-5) + (Stateminus(end-5));

RB      = StateSS(end-4) + (State(end-4));
RBminus = StateSS(end-4) + (Stateminus(end-4));

RRB      = StateSS(end-3) + (State(end-3));
RRBminus = StateSS(end-3) + (Stateminus(end-3));

% Aggregate Exogenous States
TFP   = StateSS(end-2) + (State(end-2));
TFPminus  = StateSS(end-2) + (Stateminus(end-2));

EPS_TAYLOR   = StateSS(end-1) + (State(end-1));
EPS_TAYLORminus  = StateSS(end-1) + (Stateminus(end-1));

EPS_NW       = StateSS(end) + (State(end));
EPS_NWminus  = StateSS(end) + (Stateminus(end));

%% Split the Control vector into items with names
% Controls

% Aggregate Controls (t+1)
YM = exp(Control(YMind));
N = exp(Control(Nind ));
I = exp(Control(Iind ));
C = exp(Control(Cind ));
NWE = exp(Control(NWEind ));
NWN = exp(Control(NWNind ));
MC = exp(Control(MCind ));
W  = exp(Control(Wind ));
Profit = exp(Control(Profitind ));
Y = exp(Control(Yind ));
LD = exp(Control(LDind ));
RBB  = exp(Control(RBBind ));
NU = exp(Control(NUind ));
ETA = exp(Control(ETAind ));
Z = exp(Control(Zind ));
X = exp(Control(Xind ));
PI = exp(Control(PIind));

% Aggregate Controls (t)
YMminus = exp(Controlminus(YMind));
Nminus = exp(Controlminus(Nind ));
Iminus = exp(Controlminus(Iind ));
Cminus = exp(Controlminus(Cind ));
NWEminus = exp(Controlminus(NWEind ));
NWNminus = exp(Controlminus(NWNind ));
MCminus = exp(Controlminus(MCind ));
Wminus  = exp(Controlminus(Wind ));
Profitminus = exp(Controlminus(Profitind ));
Yminus = exp(Controlminus(Yind ));
LDminus = exp(Controlminus(LDind ));
RBBminus  = exp(Controlminus(RBBind ));
NUminus = exp(Controlminus(NUind ));
ETAminus = exp(Controlminus(ETAind ));
Zminus = exp(Controlminus(Zind ));
Xminus = exp(Controlminus(Xind ));
PIminus = exp(Controlminus(PIind));

K = exp(K); Kminus = exp(Kminus);
Q = exp(Q); Qminus = exp(Qminus);
VR = exp(VR); VRminus = exp(VRminus);
NW = exp(NW); NWminus = exp(NWminus);
PHIM = exp(PHIM); PHIMminus = exp(PHIMminus);
RB = exp(RB); RBminus = exp(RBminus);
RRB = exp(RRB); RRBminus = exp(RRBminus);
TFP = exp(TFP); TFPminus  = exp(TFPminus);
EPS_TAYLOR = exp(EPS_TAYLOR); EPS_TAYLORminus = exp(EPS_TAYLORminus);
EPS_NW = exp(EPS_NW); EPS_NWminus = exp(EPS_NWminus);

%%

% f1
LHS(Kind) = VR;
RHS(Kind) = (Cminus-Nminus^(1+par.gamma)/(1+par.gamma))^(-par.xi);

% f2
LHS(Qind) = par.beta*RRB*LD;
RHS(Qind) = 1;

% f3
LHS(VRind) = LDminus;
RHS(VRind) = VR/VRminus;

% f4
LHS(NWind) = (Nminus^(par.gamma));
RHS(NWind) = (MCminus*par.alpha*YMminus/Nminus);

% f5
LHS(PHIMind) = NUminus;
RHS(PHIMind) = par.beta*LD*( (1-par.theta)*(RBB-RRB) + par.theta*X*NU );

% f6
LHS(RBind) = ETAminus;
RHS(RBind) = (1-par.theta) + par.beta*LD*( par.theta*Z*ETA );

% f7
LHS(RRBind) = PHIM;
RHS(RRBind) = ETAminus/(par.lambda-NUminus);

% f8
LHS(TFPind) = Zminus;
RHS(TFPind) = (RBBminus-RRBminus)*PHIMminus + RRBminus;

% f9
LHS(EMind) = Xminus;
RHS(EMind) = PHIM/PHIMminus*Zminus;

% f10
LHS(ENind) = K*Q;
RHS(ENind) = PHIM*NW; 

% f11
LHS(nx+YMind) = NW;
RHS(nx+YMind) = NWEminus + NWNminus;

% f12
LHS(nx+Nind) = NWEminus;
RHS(nx+Nind) = par.theta*Zminus*NWminus*EPS_NWminus;

% f13
LHS(nx+Iind) = NWNminus;
RHS(nx+Iind) = par.omega*Q*Kminus;

% f14
LHS(nx+Cind) = (RBBminus);
RHS(nx+Cind) = (MCminus*(1-par.alpha)*(YMminus./Kminus) + Q - par.delta)/Qminus;

% f15
LHS(nx+NWEind) = YMminus;
RHS(nx+NWEind) = (TFPminus*(Nminus).^(par.alpha).*Kminus.^(1-par.alpha));

% f16
LHS(nx+NWNind) = (Q);
RHS(nx+NWNind) = (par.phi*(K./Kminus-1)+1);

% f17
LHS(nx+MCind) = K - (1-par.delta)*Kminus;
RHS(nx+MCind) = Iminus;

% f18
LHS(nx+Wind) = (Yminus);
RHS(nx+Wind) = (Cminus + Iminus + par.phi/2*(K-Kminus)^2/Kminus);

% f19
LHS(nx+Profitind) = YMminus;
RHS(nx+Profitind) = Yminus + par.epsilon/2/par.kappa*Yminus*log(PIminus)^2;

% f20
LHS(nx+Yind) = log(PIminus);
RHS(nx+Yind) = par.beta*log(PI)*Y/Yminus + par.kappa*(MCminus - (par.epsilon-1)/par.epsilon);

% f21
LHS(nx+LDind) = RB;
RHS(nx+LDind) = RRB*PI;

% f22
LHS(nx+RBBind) = RB/par.RB;
RHS(nx+RBBind) = (RBminus/par.RB)^par.rho_R*(PIminus/par.PI)^((1-par.rho_R)*par.theta_pi)*EPS_TAYLORminus;

% f23
LHS(nx+NUind)  = (Profitminus);
RHS(nx+NUind) = (1-MCminus)*YMminus - YMminus.*(1/(1-par.mu))./par.kappa./2 .*log(PIminus).^2 ...
    + 1/2*par.phi*((K-Kminus).^2)./Kminus;

% f24
LHS(nx+ETAind) = Wminus;
RHS(nx+ETAind) = MCminus.* par.alpha *(YMminus./(Nminus));

% f25
LHS(nx+Zind) = log(TFP);
RHS(nx+Zind) = par.rhoS_TFP * log(TFPminus);

%f26
LHS(nx+Xind) = log(EPS_TAYLOR);
RHS(nx+Xind) = par.rhoS_MP * log(EPS_TAYLORminus);

%f27
LHS(nx+PIind) = log(EPS_NW);
RHS(nx+PIind) = par.rhoS_NW * log(EPS_NWminus);



%% Difference
% Difference = InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);(ControlSS(1:end-oc));ones(oc,1)]);
Difference = InvGamma(:,:)*((LHS-RHS)./[ones(nx,1);ones(oc,1)]);

end
