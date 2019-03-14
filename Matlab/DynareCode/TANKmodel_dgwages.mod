/*
 * Implements a basic TANK model
 */

var pi ${\pi}$ (long_name='inflation')
    y_gap ${\tilde y}$ (long_name='output gap')
    y_nat ${y^{nat}}$ (long_name='natural output')      //(in contrast to the textbook defined in deviation from steady state)
    y ${y}$ (long_name='output')
    r_real ${r^r}$ (long_name='//real interest rate')     
    i ${i}$ (long_name='nominal interrst rate')
    n ${n}$ (long_name='hours worked')
    n_R ${n_R}$ (long_name='Ricardian hours worked')
    n_K ${n_K}$ (long_name='Keynesian hours worked')
    c_R ${c_R}$ (long_name='Ricardian consumption')
    c_K ${c_K}$ (long_name='Keynesian consumption')
    w_real ${w_r}$ (long_name='//real wage')
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
    lambda ${\lambda}$ (long_name='Keynesian population share')
    Lambda ${\Lambda}$ (long_name='debt limit as multiple of steady state labor income')
    cons_share_K $\bar{C_K}$ (long_name='Keynesian consumption share')
    cons_share_R $\bar{C_R}$ (long_name='Ricardian consumption share')
    labor_share_K $\bar{N_K}$ (long_name='Keynesian labor share')
    labor_share_R $\bar{N_R}$ (long_name='Ricardian labor share')
    ;
%----------------------------------------------------------------
% Parametrization, p. 52
%----------------------------------------------------------------
set_param_value('sigma',sigma);
set_param_value('phi',phi);
set_param_value('phi_pi',phi_pi);
set_param_value('phi_y',phi_y);
set_param_value('theta',theta);
set_param_value('rho_nu',rho_nu);
set_param_value('beta',beta);
set_param_value('alpha',alpha);
set_param_value('epsilon',epsilon);
set_param_value('lambda',lambda);
set_param_value('Lambda',Lambda);

set_param_value('cons_share_K',cons_share_K);
set_param_value('cons_share_R',cons_share_R);
set_param_value('labor_share_K',labor_share_K);
set_param_value('labor_share_R',labor_share_R);

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 
//Composite parameters
#Omega=1; // We are using constant returns to scale as firms can also adjust their capital each period (1-alpha)/(1-alpha+alpha*epsilon);  //defined on page 47
#psi_n_ya=(1+phi)/(sigma*(1-alpha)+phi+alpha); //defined on page 48
#lambda_gali=(1-theta)*(1-beta*theta)/theta*Omega; //defined on page 47, named lambda there
#kappa=lambda_gali*(sigma+(phi+alpha)/(1-alpha));  //defined on page 49

//1. New Keynesian Phillips Curve eq. (21)
pi=beta*pi(+1)+kappa*y_gap;

// Change two equation below to get wages set by union 
// w_real = sigma*c_R + phi*n_R;
// w_real = sigma*c_K + phi*n_K;
w_real = sigma*y_gap + phi*n;
n_R = n_K;

c_R = c_R(+1) - 1/sigma*(i - pi(+1));
(1-Lambda*(1-beta))*c_K = w_real + n_K - Lambda*(i(-1)-pi-r_real(-1)) - beta*Lambda*r_real;

//(1-Lambda*(1-beta))*c_K = w_real + n_K  - beta*Lambda*r_real; //this version has no fisher effect

y_gap = cons_share_R*c_R + cons_share_K*c_K;
n = labor_share_R*n_R + labor_share_K*n_K;

//3. Interest Rate Rule eq. (25)
i=phi_pi*pi+phi_y*y_gap+nu;

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
stoch_simul(order = 1,irf=7) y_gap pi i r_real nu c_R c_K n_R n_K w_real;

write_latex_dynamic_model;

