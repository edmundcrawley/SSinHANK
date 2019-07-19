%% Initialize state and control vector
tic
invutil = @(u)(((1-par.xi).*u).^(1/(1-par.xi)));
invmutil = @(mu)((1./mu).^(1/par.xi));

Xss=[squeeze(sum(joint_distr,2)); ... % marginal distribution liquid
    squeeze(sum(joint_distr,1)'); ... % marginal distribution productivity
    log(par.RB);0];

par.W=par.mc;
par.PROFITS = (1-par.mc)*grid.Y;

L_new_w = L_new(:,1:mpar.nh-1); % Consider that labor supply of Es is zero

Yss=[invmutil(mutil_c(:));L_new_w(:); log(par.PI); log(grid.Y); log(par.W); 
    log(par.PROFITS);log(grid.N); log(grid.B);log(grid.C)];

%% Construct Chebyshev Polynomials to describe deviations of policy from SS
Poly=[];
Poly_L=[];
maxlevel=max([mpar.nm mpar.nh]);

% Values of Chebyshev Polynomials (1st type) at Chebyshev Nodes
Tm = cos(pi* (0:maxlevel-1)' * (linspace(0.5/mpar.nm/2,1-0.5/mpar.nm*2,mpar.nm)))';
Th = cos(pi* (0:maxlevel-1)' * (linspace(0.5/(mpar.nh-1),1-0.5/(mpar.nh-1),(mpar.nh-1))))';

for j1=1:(length(grid.h)-1)
    for j3=1:length(grid.m)
        if j1+j3<mpar.maxdim
            [TT1, TT3] = ndgrid(Tm(:,j3), [Th(:,j1); 0]);
            Poly=[Poly TT1(:).*TT3(:)];
        end
    end
end

for j2=1:length(grid.m)
    if j2<mpar.maxdim-1
        [TT1, TT3] = ndgrid(Tm(:,j2), [zeros(length(grid.h)-1,1);1]);
        Poly=[Poly TT1(:).*TT3(:)];
    end
end

for j1=1:(length(grid.h)-1)
    for j3=1:length(grid.m)
        if j1+j3<mpar.maxdim
            [TT1, TT3] = ndgrid(Tm(:,j3), [Th(:,j1)]);
            Poly_L=[Poly_L TT1(:).*TT3(:)];
        end
    end
end

InvCheb=(Poly'*Poly)\Poly';
InvCheb_L = (Poly_L'*Poly_L)\Poly_L';
%% Construct function such that perturbed marginal distributions still integrate to 1
Gamma=zeros(mpar.nm+mpar.nh,mpar.nm+mpar.nh-3);
for j=1:mpar.nm-1
    Gamma(1:mpar.nm,j)=-Xss(1:mpar.nm);
    Gamma(j,j)=1-Xss(j);
    Gamma(j,j)=Gamma(j,j) -sum(Gamma(1:mpar.nm,j));
end
bb=mpar.nm;
for j=1:mpar.nh-2
    Gamma(bb+(1:mpar.nh-1),bb+j-1)=-Xss(bb+(1:mpar.nh-1));
    Gamma(bb+j,bb-1+j)=1-Xss(bb+j);
    Gamma(bb+j,bb-1+j)=Gamma(bb+j,bb-1+j) -sum(Gamma(bb+(1:mpar.nh-1),bb-1+j));
end

%% Collect all functions used for perturbation
n1 = size(Poly); % used for controls
n1_L = size(Poly_L);
n2 = size(Gamma); %used for distributions

% Produce matrices to reduce state-space
oc = length(Yss) - (n1(1)+n1_L(1));
os = length(Xss) - (mpar.nm+mpar.nh);

InvGamma                                      = sparse((n1(1)+n1_L(1))+n2(2)+os+oc,(n1(2)+n1_L(2))+n2(2)+os+oc);

Gamma_state                                   = sparse(Gamma);
InvGamma(1:n2(1)+os,1:n2(1)+os)               = eye(n2(1)+os);

Gamma_control                                 = sparse((n1(2)+n1_L(2))+oc,(n1(2)+n1_L(2))+oc);

Gamma_control(1:n1(1),1:n1(2))                = Poly;
InvGamma(n2(2)+os+(1:n1(1)),n2(2)+os+(1:n1(2))) = InvCheb';

Gamma_control(n1(1)+(1:n1_L(1)),n1(2)+(1:n1_L(2)))  = Poly_L;
InvGamma(n2(2)+n1(1)+os+(1:n1_L(1)),n2(2)+n1(2)+os+(1:n1_L(2))) = InvCheb_L';

Gamma_control((n1(1)+n1_L(1))+(1:oc),(n1(2)+n1_L(2))+(1:oc))                  = eye(oc);
InvGamma(n2(2)+(n1(1)+n1_L(1))+os+(1:oc),n2(2)+(n1(2)+n1_L(2))+os+(1:oc))       = eye(oc);

InvGamma                                      = InvGamma';

mpar.numstates   = n2(2)+os;
mpar.numcontrols = n1(2)+n1_L(2)+oc;
State       = zeros(mpar.numstates,1);
State_m     = State;
Contr       = zeros(mpar.numcontrols,1);
Contr_m     = Contr;
