close all

IRF_distr=Gamma_state*IRF_state_sparse(1:mpar.numstates-2,1:mpar.maxlag);
% preparation
IRF_H=100*grid.h(1:end-1)*IRF_distr(mpar.nm+(1:mpar.nh-1),2:end)/par.H;
IRF_M=100*grid.m*IRF_distr((1:mpar.nm),2:end)/targets.Y;
M=100*grid.m*IRF_distr((1:mpar.nm),1:end)+grid.B;
IRF_RB=100*IRF_state_sparse(mpar.numstates-os+1,2:end);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+2,1:end-1);

Y=Output*(1+IRF_state_sparse(end-oc+2,1:end-1));
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_C=IRF_Y;
IRF_Profit=100*IRF_state_sparse(end-oc+4,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_W = IRF_N;
% IRF_B=100*IRF_state_sparse(end-3,1:end-1);
IRF_C_agg = 100*IRF_state_sparse(end-2,1:end-1);
% IRF_G=100*IRF_state_sparse(end-1,1:end-1);
% IRF_T=100*IRF_state_sparse(end,1:end-1);

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

IRF_C_by_suff = Inc_wt_MPC*IRF_Y + Earning_hetero.*IRF_Y ... 
                - Redist_elas_P/E_Y*IRF_PI/100 + Redist_elas_R/E_Y*IRF_RBREAL/100 - Hick_scaling/E_Y*IRF_RBREAL/100;

% recall that IRF_PI and IRF_RBREAL are in bp term!
%% Plotting

figure(1)
plot(1:mpar.maxlag-1,IRF_C)
hold on
plot(1:mpar.maxlag-1,IRF_C_by_suff)
ylabel('Percent', 'FontSize',12)
xlabel('Quarter', 'FontSize',12)
title('Consumption')
legend({'true IRF','IRF using sufficient stats'},'location','southeast','FontSize',12)

fprintf('M is %2.2f \n',Inc_wt_MPC)
fprintf('E(Epsilon_Y) is %2.2f \n',mean(Earning_hetero))
fprintf('Epsilon_P is %2.2f \n',mean(Redist_elas_P))
fprintf('Epsilon_R is %2.2f \n',mean(Redist_elas_R))
fprintf('Hicksian_scaling factor is %2.2f \n',mean(Hick_scaling))
%%
figure1 = figure('Position',[200 130 1200 700]);
subplot(2,3,1)
plot(1:mpar.maxlag-1,IRF_Y)
ylabel('Percent')
xlabel('Quarter')
title('Output')
subplot(2,3,2)
plot(1:mpar.maxlag-1,IRF_C)
hold on
plot(1:mpar.maxlag-1,IRF_C_agg,'--')
ylabel('Percent')
xlabel('Quarter')
title('Consumption')
subplot(2,3,3)
plot(1:mpar.maxlag-1,IRF_N)
ylabel('Percent')
xlabel('Quarter')
title('Labor')
subplot(2,3,4)
plot(1:mpar.maxlag-1,IRF_PI)
ylabel('Basis point')
xlabel('Quarter')
title('Inflation')
subplot(2,3,5)
plot(1:mpar.maxlag-1,IRF_RBREAL)
ylabel('Basis point')
xlabel('Quarter')
title('Real interest rate')
subplot(2,3,6)
plot(1:mpar.maxlag-1,IRF_S)
ylabel('Percent')
xlabel('Quarter')
title('Shock')

