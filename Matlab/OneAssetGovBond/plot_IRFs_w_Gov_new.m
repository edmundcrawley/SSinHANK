close all

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-2,1:mpar.maxlag);
% preparation
IRF_H=100*grid.h(1:end-1)*IRF_distr(mpar.nm+(1:mpar.nh-1),2:end)/par.H;
IRF_M=100*grid.m*IRF_distr((1:mpar.nm),2:end)/targets.Y;
% M=100*grid.m*IRF_distr((1:mpar.nm),1:end)+grid.B;
IRF_RB=100*IRF_state_sparse(mpar.numstates-os+1,2:end);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+2,1:end-1);

Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);

IRF_W=100*IRF_state_sparse(end-oc+3,1:end-1);
IRF_Profit=100*IRF_state_sparse(end-oc+4,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);

IRF_B=100*IRF_state_sparse(end-oc+6,1:end-1);
IRF_C_agg = 100*IRF_state_sparse(end-oc+7,1:end-1);
% IRF_G=100*IRF_state_sparse(end-oc+8,1:end-1);
IRF_X=100*IRF_state_sparse(end-oc+9,1:end-1);
IRF_MC=100*IRF_state_sparse(end-oc+10,1:end-1);

IRF_PI=100*100*IRF_state_sparse(end-oc+1,1:end-1);

PI=1+IRF_state_sparse(end-oc+1,1:end-1);
RB=par.RB+(IRF_state_sparse(mpar.numstates-os+1,2:end));
IRF_RB=100*100*(RB-par.RB);
IRF_RBREAL=100*100*(RB./PI-par.RB);


%% 2. Earning heterogeneity channel

labor_inc_share = repmat([ones(mpar.nm,mpar.nh-1),zeros(mpar.nm,1)],[1,1,length(IRF_Y)]); % nm*nh*length(IRF)
profit_inc_share = repmat([zeros(mpar.nm,mpar.nh-1),ones(mpar.nm,1)],[1,1,length(IRF_Y)]);

IRF_Yi = reshape(IRF_N+IRF_W,[1,1,length(IRF_Y)]).*labor_inc_share...
          + reshape(IRF_Profit,[1,1,length(IRF_Y)]).*profit_inc_share;
      
E_IRF_Yi = squeeze(sum(sum(repmat(joint_distr,[1,1,length(IRF_Y)]).*IRF_Yi,1),2));      


E_Y = sum(sum(joint_distr.*WW_h_mesh));
MPC_over_t = repmat(MPC_m,[1,1,length(IRF_Y)]);

Earning_hetero = sum(sum(  repmat(joint_distr,[1,1,length(IRF_Y)]).*MPC_over_t.*(IRF_Yi-reshape(IRF_Y,[1,1,length(IRF_Y)])) ...
                .*repmat(WW_h_mesh,[1,1,length(IRF_Y)])/E_Y./reshape(IRF_Y,[1,1,length(IRF_Y)])  ,1),2);
Earning_hetero = squeeze(Earning_hetero)';            

%% IRFs by sufficient statistics

% change Auclert's notations into Crawley's ones

M = Inc_wt_MPC*Output/C_agg;
Earning_hetero = Earning_hetero*Output/C_agg;
Redist_elas_P = Redist_elas_P*(1/C_agg);
Redist_elas_R = Redist_elas_R*(1/C_agg);
Hick_scaling = Hick_scaling*(1/C_agg);
            
IRF_C_by_suff = M*IRF_Y + Earning_hetero.*IRF_Y - Redist_elas_P*IRF_PI/100 ...
                + Redist_elas_R*IRF_RBREAL/100 - Hick_scaling*IRF_RBREAL/100;            

% recall that IRF_PI and IRF_RBREAL are in bp term!
%% Plotting

figure(1)
plot(1:mpar.maxlag-1,IRF_C_agg)
hold on
plot(1:mpar.maxlag-1,IRF_C_by_suff,'r--')
ylabel('Percent', 'FontSize',14)
xlabel('Quarter', 'FontSize',14)
title('Consumption')
legend({'true IRF','IRF using sufficient stats'},'location','southeast','FontSize',12)

% %% distribution
% 
% figure(3)
% plot(grid.m,sum(joint_distr,2),'LineWidth',2);xlim([0 1]);
% legend({'Distribution of wealth'},'fontsize',12)
% 
% figure(4)
% plot(grid.m,cumsum(sum(joint_distr,2)),'r--','LineWidth',2);xlim([0 1]);
% legend({'Cumulative Distribution'},'fontsize',12)

%% Compute sufficient statistics
fprintf('M is %2.2f \n',M)
fprintf('E(Epsilon_Y) is %2.2f \n',mean(Earning_hetero))
fprintf('Epsilon_P is %2.2f \n',mean(Redist_elas_P))
fprintf('Epsilon_R is %2.2f \n',mean(Redist_elas_R))
fprintf('Hicksian_scaling factor is %2.2f \n',mean(Hick_scaling))
%%
figure3 = figure('Position',[200 130 1200 700]);
subplot(3,4,1)
plot(1:mpar.maxlag-1,IRF_Y)
ylabel('Percent')
xlabel('Quarter')
title('Output')
subplot(3,4,2)
plot(1:mpar.maxlag-1,IRF_C_agg)
hold on
plot(1:mpar.maxlag-1,IRF_C_by_suff,'--')
ylabel('Percent')
xlabel('Quarter')
title('Consumption on C')
subplot(3,4,3)
plot(1:mpar.maxlag-1,IRF_N)
ylabel('Percent')
xlabel('Quarter')
title('Labor')
subplot(3,4,4)
plot(1:mpar.maxlag-1,IRF_W)
ylabel('Basis point')
xlabel('Quarter')
title('Wage')
subplot(3,4,5)
plot(1:mpar.maxlag-1,IRF_PI)
ylabel('Basis point')
xlabel('Quarter')
title('Inflation')
subplot(3,4,6)
plot(1:mpar.maxlag-1,IRF_X)
ylabel('Percent')
xlabel('Quarter')
title('Consumption on X')
subplot(3,4,7)
plot(1:mpar.maxlag-1,IRF_RBREAL)
ylabel('Basis point')
xlabel('Quarter')
title('Real interest rate')
% subplot(3,4,8)
% plot(1:mpar.maxlag-1,IRF_G)
% ylabel('Percent')
% xlabel('Quarter')
% title('Gov spending')
% subplot(3,4,9)
% plot(1:mpar.maxlag-1,IRF_T)
% ylabel('Percent')
% xlabel('Quarter')
% title('Tax')
subplot(3,4,10)
plot(1:mpar.maxlag-1,IRF_B)
ylabel('Percent')
xlabel('Quarter')
title('Bond')
subplot(3,4,11)
plot(1:mpar.maxlag-1,IRF_MC)
ylabel('Percent')
xlabel('Quarter')
title('Maginal cost')
subplot(3,4,12)
plot(1:mpar.maxlag-1,IRF_S)
ylabel('Percent')
xlabel('Quarter')
title('Shock')

