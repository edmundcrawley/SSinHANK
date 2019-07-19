close all;
clear all;clc;

tic
%% Load model and steady state values
Symbolic_SGU_model;

Symbolic_SGU_model_ss;

%%
%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

anal_deriv(f,x,y,xp,yp,2);

x0 = zeros(size(hx,1),1);      %INITIAL CONDITION for STATE VECTOR
%% Choose shock

x0(end-2) = 1; % TFP shock
% x0(end-1) = 1; % MP shock
% x0(end) = 1;   % NW shock

%% Second-order approximation
[gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);

[gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,x0);

toc


% ghx = [gx;hx];
% 
% x0 = zeros(size(hx,1),1);      %INITIAL CONDITION for STATE VECTOR
% 
% % shock = 'TFP';
% shock = 'MP';
% % shock = 'NW';
% 
% T=40;
% 
% switch(shock)
%   case('TFP')
% x0(end-2) = 1; % TFP shock
% 
% 
% IR=ir(gx,hx,x0,T);
% time=(0:size(IR,1)-1)'; 
% % Ym   L   I   C   Ne   Nn   Pm   w   profit Y Lambda   Rk  nu   eta   z   X   infl
% % K Q varrho N  phi i R a
% 
% Ym_e_a = IR(1:T,1); L_e_a = IR(1:T,2); I_e_a = IR(1:T,3); C_e_a = IR(1:T,4); 
% Ne_e_a = IR(1:T,5); Nn_e_a = IR(1:T,6); Pm_e_a = IR(1:T,7); w_e_a = IR(1:T,8); 
% profit_e_a = IR(1:T,9); Y_e_a = IR(1:T,10);Lambda_e_a = IR(1:T,11);Rk_e_a = IR(1:T,12);
% nu_e_a = IR(1:T,13); eta_e_a = IR(1:T,14); z_e_a = IR(1:T,15); X_e_a = IR(1:T,16);
% infl_e_a = IR(1:T,17); K_e_a = IR(2:T,18); Q_e_a = IR(2:T,19); varrho_e_a = IR(2:T,20);
% N_e_a = IR(2:T,21); phi_e_a = IR(2:T,22); i_e_a = IR(2:T,23); R_e_a = IR(2:T,24); ea_e_a = IR(2:T,25);
% 
% 
% figure1 = figure('Position', [50, 70, 1500, 900],'Name', 'SGU, wrt 1.0% TFP shock');
% subplot(4,5,1)
% plot(Y_e_a)
% title('Output')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,2)
% plot(C_e_a)
% title('Consumption')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,3)
% plot(I_e_a)
% title('Investment')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,4)
% plot(L_e_a)
% title('Labor supply')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,6)
% plot(profit_e_a)
% title('Profits')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,8)
% plot(K_e_a)
% title('Capital')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,9)
% plot(Pm_e_a)
% title('Marginal cost')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % subplot(4,5,9)
% % plot(B_e_a)
% % ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % title('Bond')
% subplot(4,5,10)
% plot(ea_e_a,'Linewidth',1.5)
% title('Shock')
% subplot(4,5,11)
% plot(N_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Net worth')
% subplot(4,5,12)
% plot(phi_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Leverage')
% subplot(4,5,13)
% plot(nu_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Assets')
% subplot(4,5,14)
% plot(eta_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Net Worth')
% subplot(4,5,15)
% plot(z_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Net Worth')
% subplot(4,5,16)
% plot(X_e_a)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Assets')
% subplot(4,5,17)
% plot(100*R_e_a)
% ylabel('Basis Point','Interpreter','none','FontName','arial','FontSize',10)
% title('Saving rate')
% subplot(4,5,18)
% plot(100*Rk_e_a)
% title('Borrowing rate')
% ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,19)
% plot(100*infl_e_a)
% title('Inflation')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,20)
% plot(100*Q_e_a)
% title('Capital Price')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% 
%  case('MP')
%     
% x0(end-1) = 0.1; % MP shock
% 
% 
% IR=ir(gx,hx,x0,T);
% time=(0:size(IR,1)-1)'; 
% % Ym   L   I   C   Ne   Nn   Pm   w   profit Y Lambda   Rk  nu   eta   z   X   infl
% % K Q varrho N  phi i R a
% 
% Ym_e_i = IR(1:T,1); L_e_i = IR(1:T,2); I_e_i = IR(1:T,3); C_e_i = IR(1:T,4); Ne_e_i = IR(1:T,5);
% Nn_e_i = IR(1:T,6); Pm_e_i = IR(1:T,7); w_e_i = IR(1:T,8); profit_e_i = IR(1:T,9); Y_e_i = IR(1:T,10);
% Lambda_e_i = IR(1:T,11); Rk_e_i = IR(1:T,12); nu_e_i = IR(1:T,13); eta_e_i = IR(1:T,14); z_e_i = IR(1:T,15);
% X_e_i = IR(1:T,16); infl_e_i = IR(1:T,17); K_e_i = IR(2:T,18); Q_e_i = IR(2:T,19); varrho_e_i = IR(2:T,20);
% N_e_i = IR(2:T,21); phi_e_i = IR(2:T,22); i_e_i = IR(2:T,23); R_e_i = IR(2:T,24); ei_e_i = IR(1:T,26);
% 
% 
% figure2 = figure('Position', [50, 70, 1500, 900],'Name', 'SGU, wrt 0.1% MP shock');
% subplot(4,5,1)
% plot(Y_e_i)
% title('Output')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,2)
% plot(C_e_i)
% title('Consumption')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,3)
% plot(I_e_i)
% title('Investment')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,4)
% plot(L_e_i)
% title('Labor supply')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,6)
% plot(profit_e_i)
% title('Profits')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,8)
% plot(K_e_i)
% title('Capital')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,9)
% plot(Pm_e_i)
% title('Marginal cost')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % subplot(4,5,9)
% % plot(B_e_i)
% % ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % title('Bond')
% subplot(4,5,10)
% plot(ei_e_i,'Linewidth',1.5)
% title('Shock')
% subplot(4,5,11)
% plot(N_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Net worth')
% subplot(4,5,12)
% plot(phi_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Leverage')
% subplot(4,5,13)
% plot(nu_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Assets')
% subplot(4,5,14)
% plot(eta_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Net Worth')
% subplot(4,5,15)
% plot(z_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Net Worth')
% subplot(4,5,16)
% plot(X_e_i)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Assets')
% subplot(4,5,17)
% plot(100*R_e_i)
% ylabel('Basis Point','Interpreter','none','FontName','arial','FontSize',10)
% title('Saving rate')
% subplot(4,5,18)
% plot(100*Rk_e_i)
% title('Borrowing rate')
% ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,19)
% plot(100*infl_e_i)
% title('Inflation')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,20)
% plot(100*Q_e_i)
% title('Capital Price')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% 
%    
% case('NW')
% x0(end) = -1; % TFP shock
% 
% 
% IR=ir(gx,hx,x0,T);
% time=(0:size(IR,1)-1)'; 
% % Ym   L   I   C   Ne   Nn   Pm   w   profit Y Lambda   Rk  nu   eta   z   X   infl
% % K Q varrho N  phi i R a
% 
% Ym_e_Ne = IR(1:T,1); L_e_Ne = IR(1:T,2); I_e_Ne = IR(1:T,3); C_e_Ne = IR(1:T,4); Ne_e_Ne = IR(1:T,5);
% Nn_e_Ne = IR(1:T,6); Pm_e_Ne = IR(1:T,7); w_e_Ne = IR(1:T,8); profit_e_Ne = IR(1:T,9); Y_e_Ne = IR(1:T,10);
% Lambda_e_Ne = IR(1:T,11); Rk_e_Ne = IR(1:T,12); nu_e_Ne = IR(1:T,13); eta_e_Ne = IR(1:T,14); z_e_Ne = IR(1:T,15);
% X_e_Ne = IR(1:T,16); infl_e_Ne = IR(1:T,17); K_e_Ne = IR(2:T,18); Q_e_Ne = IR(2:T,19); varrho_e_Ne = IR(2:T,20);
% N_e_Ne = IR(2:T,21); phi_e_Ne = IR(2:T,22); i_e_Ne = IR(2:T,23); R_e_Ne = IR(2:T,24); eN_e_Ne = IR(1:T,27);
% 
% 
% figure3 = figure('Position', [50, 70, 1500, 900],'Name', 'SGU, wrt -1.0% NW shock');
% subplot(4,5,1)
% plot(Y_e_Ne)
% title('Output')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,2)
% plot(C_e_Ne)
% title('Consumption')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,3)
% plot(I_e_Ne)
% title('Investment')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,4)
% plot(L_e_Ne)
% title('Labor supply')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,6)
% plot(profit_e_Ne)
% title('Profits')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,8)
% plot(K_e_Ne)
% title('Capital')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,9)
% plot(Pm_e_Ne)
% title('Marginal cost')
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % subplot(4,5,9)
% % plot(B_e_Ne)
% % ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% % title('Bond')
% subplot(4,5,10)
% plot(eN_e_Ne,'Linewidth',1.5)
% title('Shock')
% subplot(4,5,11)
% plot(N_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Net worth')
% subplot(4,5,12)
% plot(phi_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Leverage')
% subplot(4,5,13)
% plot(nu_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Assets')
% subplot(4,5,14)
% plot(eta_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Value of Net Worth')
% subplot(4,5,15)
% plot(z_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Net Worth')
% subplot(4,5,16)
% plot(X_e_Ne)
% ylabel('Percent','Interpreter','none','FontName','arial','FontSize',10)
% title('Growth of Assets')
% subplot(4,5,17)
% plot(100*R_e_Ne)
% ylabel('Basis Point','Interpreter','none','FontName','arial','FontSize',10)
% title('Saving rate')
% subplot(4,5,18)
% plot(100*Rk_e_Ne)
% title('Borrowing rate')
% ylabel('Basis points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,19)
% plot(100*infl_e_Ne)
% title('Inflation')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% subplot(4,5,20)
% plot(100*Q_e_Ne)
% title('Capital Price')
% ylabel('Basis Points','Interpreter','none','FontName','arial','FontSize',10)
% 
% end