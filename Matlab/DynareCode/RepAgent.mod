/*
 * This file implements the baseline New Keynesian model of Jordi Galí (2008): Monetary Policy, Inflation,
 * and the Business Cycle, Princeton University Press, Chapter 3
 *
 * Note that all model variables are expressed in deviations from steady state, i.e. in contrast to
 * to the chapter, both the nominal interest rate and natural output are not in log-levels, but rather mean 0
 *
 * This implementation was written by Johannes Pfeifer. In case you spot mistakes,
 * email me at jpfeifer@gmx.de
 *
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2013-15 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */

var pi ${\pi}$ (long_name='inflation')
    y_gap ${\tilde y}$ (long_name='output gap')
    y_nat ${y^{nat}}$ (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y ${y}$ (long_name='output')
    r_nat ${r^{nat}}$ (long_name='natural interest rate')
    r_real ${r^r}$ (long_name='//real interest rate')     
    i ${i}$ (long_name='nominal interrst rate')
    n ${n}$ (long_name='hours worked')
    nu ${\nu}$ (long_name='AR(1) monetary policy shock process')    
    a  ${a}$ (long_name='AR(1) technology shock process')
    ;     

varexo eps_nu ${\varepsilon_\nu}$   (long_name='monetary policy shock')
       ;

parameters alpha ${\alpha}$ (long_name='capital share')
    beta ${\beta}$ (long_name='discount factor')
    rho_nu ${\rho_{\nu}}$ (long_name='autocorrelation monetary policy shock')
    sigma ${\sigma}$ (long_name='coefficient of relative risk aversion')
    phi ${\phi}$ (long_name='Frisch elasticity')
    phi_pi ${\phi_{\pi}}$ (long_name='inflation feedback Taylor Rule')
    phi_y ${\phi_{y}}$ (long_name='output feedback Taylor Rule')
    epsilon ${\epsilon}$ (long_name='demand elasticity')
    theta ${\theta}$ (long_name='Calvo parameter')
    ;
%----------------------------------------------------------------
% Parametrization, p. 52
%----------------------------------------------------------------
sigma = 1;
phi=1;
phi_pi = 1.5;
phi_y  = .5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.99;
alpha=0.0; //1/3;
epsilon=6;

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=1; // We are using constant returns to scale as firms can also adjust their capital each period (1-alpha)/(1-alpha+alpha*epsilon);  //defined on page 47
#psi_n_ya=(1+phi)/(sigma*(1-alpha)+phi+alpha); //defined on page 48
#lambda=(1-theta)*(1-beta*theta)/theta*Omega; //defined on page 47
#kappa=lambda*(sigma+(phi+alpha)/(1-alpha));  //defined on page 49

//1. New Keynesian Phillips Curve eq. (21)
pi=beta*pi(+1)+kappa*y_gap;
//2. Dynamic IS Curve eq. (22)
y_gap=-1/sigma*(i-pi(+1)-r_nat)+y_gap(+1);
//3. Interest Rate Rule eq. (25)
i=phi_pi*pi+phi_y*y_gap+nu;
//4. Definition natural rate of interest eq. (23)
r_nat=sigma*psi_n_ya*(a(+1)-a);
//5. Definition real interest rate
r_real=i-pi(+1);
//6. Definition natural output, eq. (19)
y_nat=psi_n_ya*a;
//7. Definition output gap
y_gap=y-y_nat;
//8. Monetary policy shock
nu=rho_nu*nu(-1)+eps_nu;
//9. TFP shock
a=0;
//10. Production function (eq. 13)
y=a+(1-alpha)*n;
end;

%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------


shocks;
    var eps_nu = 0.1^2; //1 standard deviation shock of 10 basis points
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
resid(1);
steady;
check;

%----------------------------------------------------------------
% generate IRFs, replicates Figures 3.1, p. 53 (interest rate rule)
% 3.3, p. 57 (money growth rule)
%----------------------------------------------------------------
stoch_simul(order = 1,irf=7) y_gap pi i r_real nu;

write_latex_dynamic_model;

