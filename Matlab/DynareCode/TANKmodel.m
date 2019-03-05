%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'TANKmodel';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('TANKmodel.log');
M_.exo_names = 'eps_nu';
M_.exo_names_tex = '{\varepsilon_\nu}';
M_.exo_names_long = 'monetary policy shock';
M_.endo_names = 'pi';
M_.endo_names_tex = '{\pi}';
M_.endo_names_long = 'inflation';
M_.endo_names = char(M_.endo_names, 'y_gap');
M_.endo_names_tex = char(M_.endo_names_tex, '{\tilde y}');
M_.endo_names_long = char(M_.endo_names_long, 'output gap');
M_.endo_names = char(M_.endo_names, 'y_nat');
M_.endo_names_tex = char(M_.endo_names_tex, '{y^{nat}}');
M_.endo_names_long = char(M_.endo_names_long, 'natural output');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, '{y}');
M_.endo_names_long = char(M_.endo_names_long, 'output');
M_.endo_names = char(M_.endo_names, 'r_real');
M_.endo_names_tex = char(M_.endo_names_tex, '{r^r}');
M_.endo_names_long = char(M_.endo_names_long, '//real interest rate');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, '{i}');
M_.endo_names_long = char(M_.endo_names_long, 'nominal interrst rate');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, '{n}');
M_.endo_names_long = char(M_.endo_names_long, 'hours worked');
M_.endo_names = char(M_.endo_names, 'n_R');
M_.endo_names_tex = char(M_.endo_names_tex, '{n_R}');
M_.endo_names_long = char(M_.endo_names_long, 'Ricardian hours worked');
M_.endo_names = char(M_.endo_names, 'n_K');
M_.endo_names_tex = char(M_.endo_names_tex, '{n_K}');
M_.endo_names_long = char(M_.endo_names_long, 'Keynesian hours worked');
M_.endo_names = char(M_.endo_names, 'c_R');
M_.endo_names_tex = char(M_.endo_names_tex, '{c_R}');
M_.endo_names_long = char(M_.endo_names_long, 'Ricardian consumption');
M_.endo_names = char(M_.endo_names, 'c_K');
M_.endo_names_tex = char(M_.endo_names_tex, '{c_K}');
M_.endo_names_long = char(M_.endo_names_long, 'Keynesian consumption');
M_.endo_names = char(M_.endo_names, 'w_real');
M_.endo_names_tex = char(M_.endo_names_tex, '{w_r}');
M_.endo_names_long = char(M_.endo_names_long, '//real wage');
M_.endo_names = char(M_.endo_names, 'nu');
M_.endo_names_tex = char(M_.endo_names_tex, '{\nu}');
M_.endo_names_long = char(M_.endo_names_long, 'AR(1) monetary policy shock process');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, '{a}');
M_.endo_names_long = char(M_.endo_names_long, 'AR(1) technology shock process');
M_.param_names = 'alpha';
M_.param_names_tex = '{\alpha}';
M_.param_names_long = 'capital share';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, '{\beta}');
M_.param_names_long = char(M_.param_names_long, 'discount factor');
M_.param_names = char(M_.param_names, 'rho_nu');
M_.param_names_tex = char(M_.param_names_tex, '{\rho_{\nu}}');
M_.param_names_long = char(M_.param_names_long, 'autocorrelation monetary policy shock');
M_.param_names = char(M_.param_names, 'sigma');
M_.param_names_tex = char(M_.param_names_tex, '{\sigma}');
M_.param_names_long = char(M_.param_names_long, 'coefficient of relative risk aversion');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, '{\phi}');
M_.param_names_long = char(M_.param_names_long, 'Frisch elasticity');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, '{\phi_{\pi}}');
M_.param_names_long = char(M_.param_names_long, 'inflation feedback Taylor Rule');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, '{\phi_{y}}');
M_.param_names_long = char(M_.param_names_long, 'output feedback Taylor Rule');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, '{\epsilon}');
M_.param_names_long = char(M_.param_names_long, 'demand elasticity');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, '{\theta}');
M_.param_names_long = char(M_.param_names_long, 'Calvo parameter');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, '{\lambda}');
M_.param_names_long = char(M_.param_names_long, 'Keynesian population share');
M_.param_names = char(M_.param_names, 'Lambda');
M_.param_names_tex = char(M_.param_names_tex, '{\Lambda}');
M_.param_names_long = char(M_.param_names_long, 'debt limit as multiple of steady state labor income');
M_.param_names = char(M_.param_names, 'cons_share_K');
M_.param_names_tex = char(M_.param_names_tex, '\bar{C_K}');
M_.param_names_long = char(M_.param_names_long, 'Keynesian consumption share');
M_.param_names = char(M_.param_names, 'cons_share_R');
M_.param_names_tex = char(M_.param_names_tex, '\bar{C_R}');
M_.param_names_long = char(M_.param_names_long, 'Ricardian consumption share');
M_.param_names = char(M_.param_names, 'labor_share_K');
M_.param_names_tex = char(M_.param_names_tex, '\bar{N_K}');
M_.param_names_long = char(M_.param_names_long, 'Keynesian labor share');
M_.param_names = char(M_.param_names, 'labor_share_R');
M_.param_names_tex = char(M_.param_names_tex, '\bar{N_R}');
M_.param_names_long = char(M_.param_names_long, 'Ricardian labor share');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 14;
M_.param_nbr = 15;
M_.orig_endo_nbr = 14;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('TANKmodel_static');
erase_compiled_function('TANKmodel_dynamic');
M_.lead_lag_incidence = [
 0 4 18;
 0 5 0;
 0 6 0;
 0 7 0;
 1 8 0;
 2 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 19;
 0 14 0;
 0 15 0;
 3 16 0;
 0 17 0;]';
M_.nstatic = 9;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(14, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(15, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 45;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
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
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0.01;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 7;
options_.order = 1;
var_list_=[];
var_list_ = 'y_gap';
var_list_ = char(var_list_, 'pi');
var_list_ = char(var_list_, 'i');
var_list_ = char(var_list_, 'r_real');
var_list_ = char(var_list_, 'nu');
var_list_ = char(var_list_, 'c_R');
var_list_ = char(var_list_, 'c_K');
var_list_ = char(var_list_, 'n_R');
var_list_ = char(var_list_, 'n_K');
var_list_ = char(var_list_, 'w_real');
info = stoch_simul(var_list_);
save('TANKmodel_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('TANKmodel_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('TANKmodel_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('TANKmodel_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('TANKmodel_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
