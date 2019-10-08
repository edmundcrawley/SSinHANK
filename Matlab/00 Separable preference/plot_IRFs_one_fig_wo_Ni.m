close all

IRF_state = Gamma_state*IRF_state_sparse(1:mpar.numstates-os,2:end);
IRF_distr = IRF_state(1:mpar.nm,:);

IRF_control = Gamma_control*IRF_state_sparse(mpar.numstates+1:end,1:end);

IRF_C_m_h = 100*IRF_control(1:end-oc,1:end-1);

IRF_RB=100*100*IRF_state_sparse(mpar.numstates-os+1,2:end);
IRF_S=100*IRF_state_sparse(mpar.numstates-os+2,1:end-1);

Y=grid.Y*(1+IRF_state_sparse(end-oc+2,1:end-1));

IRF_PI=100*100*IRF_state_sparse(end-oc+1,1:end-1);
IRF_Y=100*IRF_state_sparse(end-oc+2,1:end-1);
IRF_W=100*IRF_state_sparse(end-oc+3,1:end-1);
IRF_Profit=100*IRF_state_sparse(end-oc+4,1:end-1);
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_B=100*IRF_state_sparse(end-oc+6,1:end-1);
IRF_C=100*IRF_state_sparse(end-oc+7,1:end-1);
IRF_G=IRF_state_sparse(end-oc+8,1:end-1);

PI=1+IRF_state_sparse(end-oc+1,1:end-1);
RB=par.RB+(IRF_state_sparse(mpar.numstates-os+1,2:end));
IRF_RBREAL=100*100*(RB./PI-par.RB);

IRF_N_m_h = -par.xi/par.gamma*IRF_C_m_h + 1/par.gamma*IRF_W.*ones(mpar.nm*mpar.nh,1);

%% 2. Earning heterogeneity channel

labor_inc_share = repmat([ones(mpar.nm,mpar.nh-1),zeros(mpar.nm,1)],[1,1,length(IRF_Y)]); % nm*nh*length(IRF)
profit_inc_share = repmat([zeros(mpar.nm,mpar.nh-1),ones(mpar.nm,1)],[1,1,length(IRF_Y)]);

IRF_Yi = (reshape(IRF_N_m_h,[mpar.nm,mpar.nh,length(IRF_Y)])+reshape(kron(IRF_W,ones(mpar.nm,mpar.nh)),[mpar.nm,mpar.nh,length(IRF_Y)])).*labor_inc_share...
          + reshape(IRF_Profit,[1,1,length(IRF_Y)]).*profit_inc_share;
      
E_IRF_Yi = squeeze(sum(sum(repmat(joint_distr,[1,1,length(IRF_Y)]).*IRF_Yi,1),2));      


E_Y = sum(sum(joint_distr.*WW_h_mesh));
MPC_over_t = repmat(MPC_m,[1,1,length(IRF_Y)]);

Earning_hetero = sum(sum(  repmat(joint_distr,[1,1,length(IRF_Y)]).*MPC_over_t.*(IRF_Yi-reshape(IRF_Y,[1,1,length(IRF_Y)])) ...
                .*repmat(WW_h_mesh,[1,1,length(IRF_Y)])/E_Y./reshape(IRF_Y,[1,1,length(IRF_Y)])  ,1),2);
Earning_hetero = squeeze(Earning_hetero)';            

%% IRFs by sufficient statistics

% change Auclert's notations into those of Crawley

M = Inc_wt_MPC*grid.Y/grid.C;
Earning_hetero = Earning_hetero*grid.Y/grid.C;
Redist_elas_P = Redist_elas_P*(1/grid.C);
Redist_elas_R = Redist_elas_R*(1/grid.C);
Hick_scaling = Hick_scaling*(1/grid.C);

IRF_AggInc = M*IRF_Y;
IRF_EarnHet = Earning_hetero.*IRF_Y;
IRF_Fisher = - Redist_elas_P*IRF_PI/100;
IRF_IRE = Redist_elas_R*IRF_RBREAL/100;
IRF_IntSubs = - Hick_scaling*IRF_RBREAL/100;

IRF_C_by_suff = IRF_AggInc + IRF_EarnHet + IRF_Fisher ...
                + IRF_IRE + IRF_IntSubs;            

%% Compute sufficient statistics
fprintf('M is %2.2f \n',M)
fprintf('E(Epsilon_Y) is %2.2f \n',mean(Earning_hetero))
fprintf('Epsilon_P is %2.2f \n',mean(Redist_elas_P))
fprintf('Epsilon_R is %2.2f \n',mean(Redist_elas_R))
fprintf('Hicksian_scaling factor is %2.2f \n',mean(Hick_scaling))
%% Plotting
% Output and Components
% figurename=['IRF_Y_theta2_' num2str(100*par.theta_pi)   '_gamma_pi_' num2str(100*par.gamma_pi) '_gamma_T_' num2str(10000*par.gamma_T) '_' casename  '_' aggrshock];

figure1=figure('Position', [500, 100, 850, 600],'Name', 'One fig');
subplot(3,3,1)
plot(1:mpar.maxlag-1,IRF_Y,'Linewidth',1.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Output')
subplot(3,3,2)
plot(1:mpar.maxlag-1,IRF_N,'Linewidth',1.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Labor supply')

subplot(3,3,3)
plot(1:mpar.maxlag-1,IRF_PI,'Linewidth',1.5)
ylabel('Basis Point','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('inflation')
subplot(3,3,4)
plot(1:mpar.maxlag-1,IRF_RB,'Linewidth',1.5);hold on;
plot(1:mpar.maxlag-1,IRF_RBREAL,'Linewidth',1.5);hold on;
ylabel('Basis Point','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
legend({'nominal','real'})
title('Saving rate')
subplot(3,3,5)
plot(1:mpar.maxlag-1,IRF_C,'Linewidth',1.5);hold on;
plot(1:mpar.maxlag-1,IRF_C_by_suff,'r--');
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Consumption (final)')
legend({'IRF C','IRF C by suff'},'location','Southeast')
subplot(3,3,6)
plot(1:mpar.maxlag-1,IRF_W,'Linewidth',1.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Wage')
subplot(3,3,7)
plot(1:mpar.maxlag-1,IRF_Profit,'Linewidth',1.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Profits')
subplot(3,3,8)
plot(1:mpar.maxlag-1,IRF_G,'Linewidth',1.5)
ylabel('\Delta from SS','Interpreter','tex','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Gov spending')
subplot(3,3,9)
plot(1:mpar.maxlag-1,IRF_S,'Linewidth',1.5)
ylabel('Percent','Interpreter','tex','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)
title('Shock')

IRF_AggInc = M*IRF_Y;
IRF_EarnHet = Earning_hetero.*IRF_Y;
IRF_Fisher = - Redist_elas_P*IRF_PI/100;
IRF_IRE = Redist_elas_R*IRF_RBREAL/100;
IRF_IntSubs = - Hick_scaling*IRF_RBREAL/100;

figure2=figure('Position', [700, 300, 450, 300],'Name', 'Decomposition');
plot(1:mpar.maxlag-1,IRF_AggInc,'y');hold on
plot(1:mpar.maxlag-1,IRF_EarnHet,'--m');hold on
plot(1:mpar.maxlag-1,IRF_Fisher,'g');hold on
plot(1:mpar.maxlag-1,IRF_IRE,'r');hold on
plot(1:mpar.maxlag-1,IRF_IntSubs,'b');hold on
plot(1:mpar.maxlag-1,IRF_C,'k','linewidth',1.5);
legend({'AggInc','EarnHet','Fisher','IRE','IntSubs','IRF C'},'location','Southeast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',12)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',12)

fprintf('IRF_C on impact: %.3f%% \n', IRF_C(1))
fprintf('IRF_C_by_suff on impact: %.3f%% \n', IRF_C_by_suff(1))

fprintf('Aggregate Income: %.1f%% \n', 100*IRF_AggInc(1)/IRF_C(1))
fprintf('Earnings Heterogeneity : %.1f%% \n', 100*IRF_EarnHet(1)/IRF_C(1))
fprintf('Fisher: %.1f%% \n', 100*IRF_Fisher(1)/IRF_C(1))
fprintf('Interest Rate Exposure: %.1f%% \n', 100*IRF_IRE(1)/IRF_C(1))
fprintf('Intertemporal Substitution: %.1f%% \n', 100*IRF_IntSubs(1)/IRF_C(1))
fprintf('Error: %.1f%% \n', 100*(IRF_C(1)-IRF_C_by_suff(1))/IRF_C(1))


