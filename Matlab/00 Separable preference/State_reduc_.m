tic
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

marginal_hminus=squeeze(sum(joint_distr,1)');

Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
     squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
     log(par.phiB); log(par.Q); log(par.RB);log(grid.K);log(par.NetWorth);0]; %  0 is the SS of aggr exogenous shock

par.CB = (1-par.theta)*par.z*par.NetWorth-par.Nn;
par.C = C_agg;
par.X = c_guess(:)'*joint_distr(:);

Yss=[invmutil(mutil_c(:)); log(par.Q);log(par.PI); log(Output);...
    log(par.W); log(par.RBB-1);log(par.PROFITS); log(grid.N);...
    log(grid.B);log(par.nuB);log(par.etaB);log(C_agg);log(par.mu);...
    log(par.Nn);log(par.z);log(par.x);log(par.delta*grid.K);log(par.Ne);...
    log(par.CB);log(par.RB);log(grid.D);log(grid.ll);log(par.X);log(par.lambda)];


[Poly,InvCheb,Gamma] = createSparseBasis(grid,mpar,mpar.maxdim,Xss);

% Collect all functions used for perturbation
n1 = size(Poly); % used for controls
n2 = size(Gamma); %used for distributions

% Produce matrices to reduce state-space
oc = length(Yss) - n1(1);
os = length(Xss) - (mpar.nm);

InvGamma                                      = sparse(1*n1(1)+n2(2)+os+oc,1*n1(2)+n2(2)+os+oc);

Gamma_state                                   = sparse(Gamma);
InvGamma(1:n2(1)+os,1:n2(1)+os)                 = eye(n2(1)+os);

Gamma_control                                 = sparse(1*n1(1)+oc,1*n1(2)+oc);
Gamma_control(1:n1(1),1:n1(2))                = Poly;
InvGamma(n2(2)+os+(1:n1(1)),n2(2)+os+(1:n1(2))) = InvCheb';

Gamma_control(1*n1(1)+(1:oc),1*n1(2)+(1:oc))                  = eye(oc);
InvGamma(n2(2)+1*n1(1)+os+(1:oc),n2(2)+1*n1(2)+os+(1:oc))       = eye(oc);

InvGamma                                      = InvGamma';

mpar.numstates   = n2(2)+os;
mpar.numcontrols = n1(2)+oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
disp('Computing system for SS.');
toc
