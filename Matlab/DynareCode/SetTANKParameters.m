%% Sets parameters for the paper

sigma = 1.0;
phi=1.0;
phi_pi = 1.5;
phi_y  = 0.0;%.5/4;
theta=2/3;
rho_nu =0.0;
beta = 0.97;
alpha=0.33;
epsilon=6;
lambda = 0.2;
Lambda = 0.0;
labor_share_K = lambda;
labor_share_R = 1-lambda;
% parameters for capital model
delta = 0.1;
psi_c = 0;

% Make table of calibration
table_dir = '.\Tables\';
calib_table = fopen([table_dir,'calib_table.tex'],'wt');
fprintf(calib_table, '  \\begin{table}\n');
fprintf(calib_table, '\\begin{center}\n');
fprintf(calib_table, '    \\caption{Baseline Calibration}\\label{table:calibration}\n');
fprintf(calib_table, '\\begin{tabular}{lcl}  \n');
fprintf(calib_table, '\\\\ \\toprule  \n');
fprintf(calib_table, '$\\sigma$ & %1.1f & Inverse EIS \\\\ \n', sigma);
fprintf(calib_table, '$\\psi$ & %1.1f & Inverse Frisch Elasticity \\\\ \n', phi);
fprintf(calib_table, '$\\phi_{\\pi}$ & %1.1f & Taylor Rule Coefficient \\\\ \n', phi_pi);
fprintf(calib_table, '$\\theta$ & %1.3f & Calvo stickiness parameter \\\\ \n', theta);
fprintf(calib_table, '$\\beta$ & %1.1f &  Discount Factor\\\\ \n', beta);
fprintf(calib_table, '$\\alpha$ & %1.2f &  Capital Share \\\\ \n', alpha);
fprintf(calib_table, '$\\varepsilon$ & %1.1f &  Elasticity of sub. between goods\\\\ \n', epsilon);
fprintf(calib_table, '$\\lambda$ & %1.1f & Share of Keynesian Households \\\\ \n', lambda);
fprintf(calib_table, '$\\Omega$ & %1.1f & Keynesian Debt as Share of Income \\\\ \n', Lambda);
fprintf(calib_table, '$\\delta$ & %1.1f & Depreciation (capital model only) \\\\ \n', delta);
fprintf(calib_table, '\\\\ \\bottomrule \n ');
fprintf(calib_table, '\\end{tabular}\n');
fprintf(calib_table, '\\end{center}\n');
fprintf(calib_table, '\\end{table}\n');
fclose(calib_table);
