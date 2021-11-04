%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'TANKmodel_dgwages';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('TANKmodel_dgwages.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'eps_nu'};
M_.exo_names_tex(1) = {'{\varepsilon_\nu}'};
M_.exo_names_long(1) = {'monetary policy shock'};
M_.endo_names = cell(14,1);
M_.endo_names_tex = cell(14,1);
M_.endo_names_long = cell(14,1);
M_.endo_names(1) = {'pi'};
M_.endo_names_tex(1) = {'{\pi}'};
M_.endo_names_long(1) = {'inflation'};
M_.endo_names(2) = {'y_gap'};
M_.endo_names_tex(2) = {'{\tilde y}'};
M_.endo_names_long(2) = {'output gap'};
M_.endo_names(3) = {'y_nat'};
M_.endo_names_tex(3) = {'{y^{nat}}'};
M_.endo_names_long(3) = {'natural output'};
M_.endo_names(4) = {'y'};
M_.endo_names_tex(4) = {'{y}'};
M_.endo_names_long(4) = {'output'};
M_.endo_names(5) = {'r_real'};
M_.endo_names_tex(5) = {'{r^r}'};
M_.endo_names_long(5) = {'//real interest rate'};
M_.endo_names(6) = {'i'};
M_.endo_names_tex(6) = {'{i}'};
M_.endo_names_long(6) = {'nominal interrst rate'};
M_.endo_names(7) = {'n'};
M_.endo_names_tex(7) = {'{n}'};
M_.endo_names_long(7) = {'hours worked'};
M_.endo_names(8) = {'n_R'};
M_.endo_names_tex(8) = {'{n_R}'};
M_.endo_names_long(8) = {'Ricardian hours worked'};
M_.endo_names(9) = {'n_K'};
M_.endo_names_tex(9) = {'{n_K}'};
M_.endo_names_long(9) = {'Keynesian hours worked'};
M_.endo_names(10) = {'c_R'};
M_.endo_names_tex(10) = {'{c_R}'};
M_.endo_names_long(10) = {'Ricardian consumption'};
M_.endo_names(11) = {'c_K'};
M_.endo_names_tex(11) = {'{c_K}'};
M_.endo_names_long(11) = {'Keynesian consumption'};
M_.endo_names(12) = {'w_real'};
M_.endo_names_tex(12) = {'{w_r}'};
M_.endo_names_long(12) = {'//real wage'};
M_.endo_names(13) = {'nu'};
M_.endo_names_tex(13) = {'{\nu}'};
M_.endo_names_long(13) = {'AR(1) monetary policy shock process'};
M_.endo_names(14) = {'a'};
M_.endo_names_tex(14) = {'{a}'};
M_.endo_names_long(14) = {'AR(1) technology shock process'};
M_.endo_partitions = struct();
M_.param_names = cell(15,1);
M_.param_names_tex = cell(15,1);
M_.param_names_long = cell(15,1);
M_.param_names(1) = {'alpha'};
M_.param_names_tex(1) = {'{\alpha}'};
M_.param_names_long(1) = {'capital share'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'discount factor'};
M_.param_names(3) = {'rho_nu'};
M_.param_names_tex(3) = {'{\rho_{\nu}}'};
M_.param_names_long(3) = {'autocorrelation monetary policy shock'};
M_.param_names(4) = {'sigma'};
M_.param_names_tex(4) = {'{\sigma}'};
M_.param_names_long(4) = {'coefficient of relative risk aversion'};
M_.param_names(5) = {'phi'};
M_.param_names_tex(5) = {'{\phi}'};
M_.param_names_long(5) = {'Frisch elasticity'};
M_.param_names(6) = {'phi_pi'};
M_.param_names_tex(6) = {'{\phi_{\pi}}'};
M_.param_names_long(6) = {'inflation feedback Taylor Rule'};
M_.param_names(7) = {'phi_y'};
M_.param_names_tex(7) = {'{\phi_{y}}'};
M_.param_names_long(7) = {'output feedback Taylor Rule'};
M_.param_names(8) = {'epsilon'};
M_.param_names_tex(8) = {'{\epsilon}'};
M_.param_names_long(8) = {'demand elasticity'};
M_.param_names(9) = {'theta'};
M_.param_names_tex(9) = {'{\theta}'};
M_.param_names_long(9) = {'Calvo parameter'};
M_.param_names(10) = {'lambda'};
M_.param_names_tex(10) = {'{\lambda}'};
M_.param_names_long(10) = {'Keynesian population share'};
M_.param_names(11) = {'Lambda'};
M_.param_names_tex(11) = {'{\Lambda}'};
M_.param_names_long(11) = {'debt limit as multiple of steady state labor income'};
M_.param_names(12) = {'cons_share_K'};
M_.param_names_tex(12) = {'\bar{C_K}'};
M_.param_names_long(12) = {'Keynesian consumption share'};
M_.param_names(13) = {'cons_share_R'};
M_.param_names_tex(13) = {'\bar{C_R}'};
M_.param_names_long(13) = {'Ricardian consumption share'};
M_.param_names(14) = {'labor_share_K'};
M_.param_names_tex(14) = {'\bar{N_K}'};
M_.param_names_long(14) = {'Keynesian labor share'};
M_.param_names(15) = {'labor_share_R'};
M_.param_names_tex(15) = {'\bar{N_R}'};
M_.param_names_long(15) = {'Ricardian labor share'};
M_.param_partitions = struct();
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
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 14;
M_.eq_nbr = 14;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
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
M_.dynamic_tmp_nbr = [1; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'pi' ;
  2 , 'name' , 'w_real' ;
  3 , 'name' , 'n_R' ;
  4 , 'name' , 'c_R' ;
  5 , 'name' , '5' ;
  6 , 'name' , 'y_gap' ;
  7 , 'name' , 'n' ;
  8 , 'name' , 'i' ;
  9 , 'name' , 'r_real' ;
  10 , 'name' , 'y_nat' ;
  11 , 'name' , '11' ;
  12 , 'name' , 'nu' ;
  13 , 'name' , 'a' ;
  14 , 'name' , 'y' ;
};
M_.mapping.pi.eqidx = [1 4 5 8 9 ];
M_.mapping.y_gap.eqidx = [1 2 6 8 11 ];
M_.mapping.y_nat.eqidx = [10 ];
M_.mapping.y.eqidx = [11 14 ];
M_.mapping.r_real.eqidx = [5 9 ];
M_.mapping.i.eqidx = [4 5 8 9 ];
M_.mapping.n.eqidx = [2 7 14 ];
M_.mapping.n_R.eqidx = [3 7 ];
M_.mapping.n_K.eqidx = [3 5 7 ];
M_.mapping.c_R.eqidx = [4 6 ];
M_.mapping.c_K.eqidx = [5 6 ];
M_.mapping.w_real.eqidx = [2 5 ];
M_.mapping.nu.eqidx = [8 12 ];
M_.mapping.a.eqidx = [13 ];
M_.mapping.eps_nu.eqidx = [12 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [5 6 13 ];
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
M_.endo_trends = struct('deflator', cell(14, 1), 'log_deflator', cell(14, 1), 'growth_factor', cell(14, 1), 'log_growth_factor', cell(14, 1));
M_.NNZDerivatives = [41; 0; -1; ];
M_.static_tmp_nbr = [1; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
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
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0.01;
resid(1);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 7;
options_.order = 1;
var_list_ = {'y_gap';'pi';'i';'r_real';'nu';'c_R';'c_K';'n_R';'n_K';'w_real'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('TANKmodel_dgwages_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('TANKmodel_dgwages_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
