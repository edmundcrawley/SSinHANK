%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for TANK model in Dynare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
addpath(genpath('c:\dynare'));
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Basic TANK model, experiment with different parameters
run 'SetTANKParameters.m';
cons_share_K = (1-Lambda*(1-beta))*lambda*(1-alpha)*(epsilon-1)/epsilon;
cons_share_R = 1-cons_share_K;
% Loop over parameters
dynare TANKmodel_dgwages.mod noclearall;
options_.nograph = 1;

%% run different experiments
figure_dir = '.\Figures\';
param_to_vary = [string('lambda'),string('Lambda')];
title_string = [string('Transmission Channels to a 1% Interest Rate Decline'),string('Transmission Channels to a 1% Interest Rate Decline'),string('Transmission Channels to a 1% Interest Rate Decline');string('\sigma = 1'),string('\sigma = 2'),string('\sigma = 3')];
save_name = [string('ProportionKeynesian'),string('KeynesianDebt')];
plot_xlabel = [string('Proportion of Keynesian Households'),string('Keynesian Household Debt Level')];
param_loop.(char(param_to_vary(1))) = 0.0:0.01:0.3;
param_loop.(char(param_to_vary(2))) = 0.0:0.05:2;
for sigma = 1:3
    set_param_value('sigma',sigma);
    for j=1:length(param_to_vary)
        set_param_value('lambda',0.2);
        set_param_value('Lambda',0.0);
        param_array = param_loop.(char(param_to_vary(j)));
        Transmission_Channels = zeros(length(param_loop),7);
        suff_stats = zeros(length(param_array),4);
        YRP_changes = zeros(length(param_array),3);
        checks = zeros(length(param_array),3);
        for i=1:length(param_array)
            set_param_value(char(param_to_vary(j)),param_array(i));
            run 'CalcTransmissionChannels_TANKmodel.m';
            Transmission_Channels(i,1) = agg_inc*nominal_i_scale;
            Transmission_Channels(i,2) = het_inc*nominal_i_scale;
            Transmission_Channels(i,3) = ire*nominal_i_scale;
            Transmission_Channels(i,4) = fisher*nominal_i_scale;
            Transmission_Channels(i,5) = ies*nominal_i_scale;
            Transmission_Channels(i,6) = dC_C_TANK*nominal_i_scale;
            Transmission_Channels(i,7) = param_array(i);
            suff_stats(i,1) = Inc_wt_MPC_TANK;
            suff_stats(i,2) = Elas_R_TANK;
            suff_stats(i,3) = Elas_P_TANK;
            suff_stats(i,4) = Elas_EIS_TANK;
            YRP_changes(i,1) = dY_Y_TANK*nominal_i_scale;
            YRP_changes(i,2) = dR_R_TANK*nominal_i_scale;
            YRP_changes(i,3) = dP_P_TANK*nominal_i_scale;
            checks(i,1) = error_TANK;
            checks(i,2) = check_R;
            checks(i,3) = check_K;
        end
        cum_trans_channels = cumsum(Transmission_Channels(:,1:5),2);
        figure;
        area(param_array,Transmission_Channels(:,[5,3,4,1,2]));
        colormap(linspecer);
        legend('Intertemporal Substitution', 'Interest Rate Exposure', 'Fisher','Aggregate Income', 'Heterogenous Income','location','NorthWest');
        legend('boxoff');
        xlabel(char(plot_xlabel(j)));
        ylabel('Consumption Change');
        title(title_string(j,sigma));
        axis([ param_array(1) param_array(end) 0 3])
        if j==1
            hold on
            plot([0.2 0.2],ylim,'Color','k','LineStyle','--');
        end
        saveas(gcf,[figure_dir,char(save_name(j)),'_sigma',num2str(sigma),'.eps'],'epsc');
    end
end



%%
% Now the TANK_capital model
run 'SetTANKParameters.m';
% Calc steady state share of labor and consumption of each type
cons_share_to_labor_share_K = (1-Lambda*(1-beta))*(epsilon-1)/epsilon*(1-alpha);
invest_share = delta*alpha*(epsilon-1)/(epsilon*(1/beta - (1-delta)));
cons_share_K = lambda*cons_share_to_labor_share_K*(1-invest_share);
cons_share_R = 1-invest_share-cons_share_K;
dynare 'TANK_capital_model.mod' noclearall;
options_.nograph = 1;

%% run different experiments
figure_dir = '.\Figures\';
psi_loop = [0,1,3,10000];
Transmission_Channels = zeros(length(psi_loop),7);
suff_stats = zeros(length(psi_loop),4);
YRP_changes = zeros(length(psi_loop),3);
IRF_i = zeros(length(psi_loop),6);
IRF_c_K = zeros(length(psi_loop),6);
IRF_c_R = zeros(length(psi_loop),6);
IRF_k = zeros(length(psi_loop),6);
IRF_q = zeros(length(psi_loop),6);
IRF_r_real = zeros(length(psi_loop),6);
checks = zeros(length(psi_loop),3);
for i=1:length(psi_loop)
    set_param_value('psi_c',psi_loop(i));
    run 'CalcTransmissionChannels_TANK_capital_model.m';
    Transmission_Channels(i,1) = agg_inc*nominal_i_scale;
    Transmission_Channels(i,2) = het_inc*nominal_i_scale;
    Transmission_Channels(i,3) = ire*nominal_i_scale;
    Transmission_Channels(i,4) = fisher*nominal_i_scale;
    Transmission_Channels(i,5) = ies*nominal_i_scale;
    Transmission_Channels(i,6) = dC_C_TANK_Capital*nominal_i_scale;
    Transmission_Channels(i,7) = psi_loop(i);
    suff_stats(i,1) = Inc_wt_MPC_TANK_Capital;
    suff_stats(i,2) = Elas_R_TANK_Capital;
    suff_stats(i,3) = Elas_P_TANK_Capital;
    suff_stats(i,4) = Elas_EIS_TANK_Capital;
    YRP_changes(i,1) = dY_Y_TANK_Capital*nominal_i_scale;
    YRP_changes(i,2) = dR_R_TANK_Capital*nominal_i_scale;
    YRP_changes(i,3) = dP_P_TANK_Capital*nominal_i_scale;
    checks(i,1) = 100*error_TANK_Capital/dC_C_TANK_Capital;
    checks(i,2) = 100*check_R/dC_C_TANK_Capital_R;
    checks(i,3) = 100*check_K/dC_C_TANK_Capital_K;
    IRF_i(i,:) = TANK_Capital_irfs.i_eps_nu(1:6)*nominal_i_scale;
    IRF_c_R(i,:) = TANK_Capital_irfs.c_R_eps_nu(1:6)*nominal_i_scale;
    IRF_c_K(i,:) = TANK_Capital_irfs.c_K_eps_nu(1:6)*nominal_i_scale;
    IRF_k(i,:) = TANK_Capital_irfs.k_eps_nu(1:6)*nominal_i_scale;
    IRF_q(i,:) = TANK_Capital_irfs.q_eps_nu(1:6)*nominal_i_scale;
    IRF_r_real(i,:) = TANK_Capital_irfs.r_real_eps_nu(1:6)*nominal_i_scale;
end
%Plot IRF for nominal interest rate
figure;
plot(IRF_i','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','SouthEast');
legend('boxoff');
xlabel('Time (years)');
title('Nominal Interest Rate Path');
ylabel('Nominal Interest Rate');
saveas(gcf,[figure_dir,'TANK_capital_IRF_i.eps'],'epsc');

%Plot IRF for c_K
figure;
plot(IRF_c_K','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','NorthEast');
legend('boxoff');
xlabel('Time (years)');
title('Keynesian Consumption Path');
ylabel('Consumption');
saveas(gcf,[figure_dir,'TANK_capital_IRF_c_K.eps'],'epsc');

%Plot IRF for c_R
figure;
plot(IRF_c_R','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','NorthEast');
legend('boxoff');
xlabel('Time (years)');
title('Ricardian Consumption Path');
ylabel('Consumption');
saveas(gcf,[figure_dir,'TANK_capital_IRF_c_R.eps'],'epsc');

%Plot IRF for Capital
figure;
plot(IRF_k','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','NorthEast');
legend('boxoff');
xlabel('Time (years)');
title('Capital Path');
ylabel('Capital');
saveas(gcf,[figure_dir,'TANK_capital_IRF_k.eps'],'epsc');

%Plot IRF for Tobin's q
figure;
plot(IRF_q','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','NorthEast');
legend('boxoff');
xlabel('Time (years)');
title('Tobin q Path');
ylabel('Tobin q');
saveas(gcf,[figure_dir,'TANK_capital_IRF_q.eps'],'epsc');

%Plot IRF for Real Interest Rate
figure;
plot(IRF_r_real','LineWidth',1.5);
xticks(1:6)
colormap(linspecer);
legend('\psi_c = 0','\psi_c = 1','\psi_c = 3','\psi_c = \infty','location','NorthEast');
legend('boxoff');
xlabel('Time (years)');
title('Real Interest Rate Path');
ylabel('Real Interest Rate');
saveas(gcf,[figure_dir,'TANK_capital_IRF_r_real.eps'],'epsc');

%Make table of Auclert statistic errors

% Make table of calibration
table_dir = '.\Tables\';
error_table = fopen([table_dir,'error_table.tex'],'wt');
fprintf(error_table, '  \\begin{table}\n');
fprintf(error_table, '\\begin{center}\n');
fprintf(error_table, '    \\caption{Percentage Error of Decomposition}\\label{table:error}\n');
fprintf(error_table, '\\begin{tabular}{cccc}  \n');
fprintf(error_table, '$\\psi_c$ & Total Consumption & Ricardian Consumption & Keynesian Consumption \n');
fprintf(error_table, '\\\\ \\toprule  \n');
fprintf(error_table, '%1.0f & %1.1f \\%%  & %1.1f \\%% & 0.0 \\%% \\\\ \n', psi_loop(1), round(checks(1,1:2),10));
fprintf(error_table, '%1.0f & %1.1f \\%%  & %1.1f \\%% & 0.0 \\%% \\\\ \n', psi_loop(2), round(checks(2,1:2),10));
fprintf(error_table, '%1.0f & %1.1f \\%%  & %1.1f \\%% & 0.0 \\%% \\\\ \n', psi_loop(3), round(checks(3,1:2),10));
fprintf(error_table, '$\\infty$ & %1.1f \\%%  & %1.1f \\%% & 0.0 \\%% \\\\ \n', round(checks(4,1:2),10));
fprintf(error_table, '\\\\ \\bottomrule \n ');
fprintf(error_table, '\\end{tabular}\n');
fprintf(error_table, '\\end{center}\n');
fprintf(error_table, '\\end{table}\n');
fclose(error_table);




