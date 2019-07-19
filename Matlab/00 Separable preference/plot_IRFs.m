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
IRF_N=100*IRF_state_sparse(end-oc+5,1:end-1);
IRF_PI=100*100*IRF_state_sparse(end-oc+1,1:end-1);

PI=1+IRF_state_sparse(end-oc+1,1:end-1);
RB=par.RB+(IRF_state_sparse(mpar.numstates-os+1,2:end));
IRF_RB=100*100*(RB-par.RB);
IRF_RBREAL=100*100*(RB./PI-par.RB);

%% Plotting
% Output and Components
figurename=['IRF_Y_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_Y,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black')
set(gca, 'FontName','arial','FontSize',25);
set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_C_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_C,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])

%%

figurename=['IRF_M_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_M,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
ylim([-1 1])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_H_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_H,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
ylim([-1 1])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_S_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_S,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])
%%
figurename=['IRF_RBPI_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_RB,'--','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_RBREAL,'LineWidth',3.5)
hold on
legend({'nominal','real'},'Location','NorthEast')
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_RB_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_RB,'LineWidth',3.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_PI_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_PI,'LineWidth',3.5)
ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])

figurename=['IRF_N_thetapi_' num2str(100*par.theta_pi)   '_rhoR_' num2str(100*par.rho_R) '_' casename  '_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 750 500])
plot(1:mpar.maxlag-1,IRF_N,'LineWidth',3.5)
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',25)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',25)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black') 
set(gca, 'FontName','arial','FontSize',25); set(gca,'Position',[0.2 0.2 0.75 0.75])
printpdf(gcf,['latex/' figurename])
%%
