%NEOCLASSICAL_MODEL_SS.M
% function [priceadj, betta, sig, varphi, theta, alfa, eta_i, epsilon, kappa_pi, rho_i, rho_a, sigma_a, kappa, omega, lambda, chi, delta_c, i_ss,...
%           Ym, L, I, C, Ne, Nn, Pm, w, profit,  K, Q, varrho, N,  phi, i, R, a,  Y, Lambda, Rk, nu, eta, z, x,infl]=GK_model_ss
%This program produces the the deep structural parameters and computes the steady state of the simple neoclassical growth model described in section 2.1 of ``Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function,'' by Stephanie Schmitt-Grohe and Martin Uribe, (2001). 
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001, revised 22-Oct-2004

theta=0.972;
lambda=0.3815;
omega=0.002;
betta=0.99000000;
sig=2.00000000;
varphi=1;
alfa=0.33000000;
eta_i=11.4;
epsilon=20;
kappa_pi=1.25000000;
rho_i=0.80000000;
rho_e_i = 0;
rho_e_N = 0;
rho_a=0.95000000;
sigma_a=0.01000000;
priceadj = 3/4; 
kappa=(1-priceadj)*(1-priceadj*betta)/priceadj;
chi=1;
delta_c=0.02;
i_ss = 1/betta;


Ym = 1.77758; Ymp=Ym;
L  = 0.662903; Lp=L;
I  =  0.12868; Ip = I;
C  =  1.56406; Cp = C;
Ne = 2.62282; Nep= Ne;
Nn = -2.17391; Nnp = Nn;
Pm = -0.0512933; Pmp = Pm;
w  =  0.662904; wp = w;
profit = -1.21814; profitp = profit;
K =  4.0407; Kp = K;
Q =  0; Qp = Q;
varrho =  -2.12637; varrhop = varrho;
N =  2.63104; Np = N;
phi =1.40966; phip = phi;
i =  0.0100503; ip = i;
R =  0.0100503; Rp = R;
aa =  0; aap = aa;
Y =  1.77758; Yp = Y;
Lambda = 0; Lambdap = Lambda;
Rk = 0.0125329; Rkp = Rk;
nu = -5.56132; nup = nu;
eta =0.435892; etap = eta;
z =  0.0201766; zp = z;
X =  0.0201766; Xp = X;
infl = 5.13567e-09; inflp = infl;
e_i = 0; e_ip = e_i;
e_N = 0; e_Np = e_N;
e_a = 0; e_ap = e_a;
% main_RANK;
% 
% sig=par.xi;
% varphi = par.gamma;
% varrho=((Output-par.delta*grid.K)-grid.N^(1+varphi)/(1+varphi))^(-sig);
% 
% Ym = log(Output); Ymp=Ym;
% L  = log(grid.N); Lp=L;
% I  = log(par.delta*grid.K); Ip = I;
% C  = log(Ym-par.delta*grid.K); Cp = C;
% Ne = log(par.Ne); Nep= Ne;
% Nn = log(par.Nn); Nnp = Nn;
% Pm = log(par.mu); Pmp = Pm;
% w  = log(par.W); wp = w;
% profit = log(par.PROFITS); profitp = profit;
% K =  log(grid.K); Kp = K;
% Q =  0; Qp = Q;
% varrho = log(-2.12637); varrhop = varrho;
% N =  log(par.NetWorth); Np = N;
% phi = log(par.phiB); phip = phi;
% i =  log(par.RB); ip = i;
% R =  log(par.RB); Rp = R;
% a =  0; ap = a;
% Y =  Ym; Yp = Y;
% Lambda = 0; Lambdap = Lambda;
% Rk = log(par.Rk); Rkp = Rk;
% nu = log(par.nuB); nup = nu;
% eta = log(par.etaB); etap = eta;
% z = log(par.z); zp = z;
% X = log(par.x); Xp = X;
% infl = 0; inflp = infl;



% end