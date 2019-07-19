%Define parameters
syms priceadj betta sig varphi theta alfa eta_i epsilon kappa_pi rho_i rho_e_i rho_e_N rho_a sigma_a kappa omega lambda chi delta_c i_ss

%Define variables 
% syms Ym   L   I   C   Ne   Nn   Pm   w   profit    K   Q   varrho   N    phi   i   R   aa  e_a Y   Lambda   Rk  nu   eta   z   X   infl e_i e_N...
%      Ymp  Lp  Ip  Cp  Nep  Nnp  Pmp  wp  profitp   Kp  Qp  varrhop  Np   phip  ip  Rp  aap e_ap Yp  Lambdap  Rkp nup  etap  zp  Xp  inflp e_ip e_Np
 
syms Ym   L   I   C   Ne   Nn   Pm   w   profit    K   Q   varrho   N    phi   i   R   aa   Y   Lambda   Rk  nu   eta   z   X   infl e_i e_N...
     Ymp  Lp  Ip  Cp  Nep  Nnp  Pmp  wp  profitp   Kp  Qp  varrhop  Np   phip  ip  Rp  aap  Yp  Lambdap  Rkp nup  etap  zp  Xp  inflp e_ip e_Np 
 
%Write equations fi, i=1:3
f1 = varrhop - (C-L^(1+varphi)/(1+varphi))^(-sig);
f2 = betta*Rp*Lambdap  -  1;
f3 = Lambda - varrhop/varrho;
f4 = chi*L^varphi - Pm*(1-alfa)*Ym/L;
f5 = nu - (1-theta)*betta*Lambdap*(Rkp - Rp) -  betta*Lambdap*theta*Xp*nup;
f6 = eta - (1-theta) - betta*Lambdap*theta*zp*etap;
f7 = phip - eta/(lambda-nu);
f8 = z - (Rk-R)*phi - R;
f9 = X - (phip/phi)*z;
f10 = Qp*Kp - phip*Np;
f11 = Np - Ne - Nn;
f12 = Ne - theta*z*N*e_N;
f13 = Nn - omega*Qp*K;
f14 = Rk - (Pm*alfa*Ym/K + (Qp-delta_c))/Q;
f15 = Ym - aa*(K)^alfa*L^(1-alfa);
f16 = Qp - 1 - eta_i*(Kp-K)/K;
f17 = Kp - (1-delta_c)*K - I;
f18 = Y - C - I - eta_i/2*((Kp-K))^2/K;
f19 = Ym - Y - epsilon/(2*kappa)*Y*log(infl)^2;
f20 = log(infl) - betta*log(inflp)*Yp/Y - kappa*(Pm - (epsilon-1)/(epsilon));
f21 = ip - Rp*inflp;
f22 = ip/i_ss - (i/i_ss)^rho_i*(infl)^(kappa_pi*(1-rho_i))*e_i;
% f22 = log(ip) - log(i_ss) - rho_i*(log(i)-log(i_ss)) - kappa_pi*(1-rho_a)*log(infl);
f23 = profit - (1-Pm)*Ym + epsilon/2/kappa*Ym*log(infl)^2 - eta_i/2*(Kp-K)^2/K;
f24 = w - Pm*(1-alfa)*Ym/L;
% f25 = log(aap) - rho_a * log(aa) - e_a;
f25 = log(aap) - rho_a * log(aa);
f26 = log(e_ip) - rho_e_i*log(e_i);
f27 = log(e_Np) - rho_e_N*log(e_N);
% f28 = e_ap;
%Create function f
f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;f19;f20;...
     f21;f22;f23;f24;f25;f26;f27];

% Define the vector of controls, y, and states, x
x = [K Q varrho N  phi i R aa e_i e_N];
y = [Ym   L   I   C   Ne   Nn   Pm   w   profit Y Lambda   Rk  nu   eta   z   X   infl];
xp = [Kp  Qp  varrhop  Np   phip  ip  Rp  aap e_ip e_Np];
yp = [Ymp  Lp  Ip  Cp  Nep  Nnp  Pmp  wp  profitp Yp  Lambdap  Rkp nup  etap  zp  Xp  inflp];

% x = [K Q varrho phi i R a e_i e_N];
% y = [N  Ym   L   I   C   Ne   Nn   Pm   w   profit Y Lambda   Rk  nu   eta   z   X   infl];
% xp = [Kp  Qp  varrhop  phip  ip  Rp  ap e_ip e_Np];
% yp = [Np Ymp  Lp  Ip  Cp  Nep  Nnp  Pmp  wp  profitp Yp  Lambdap  Rkp nup  etap  zp  Xp  inflp];

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));
%if line 36 gives an error (whichh it will for some versions of Matlab), percentage line 36 out and instead use line 38.
%f = subs(f, [x,y,xp,yp], transpose(exp([x,y,xp,yp])));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);