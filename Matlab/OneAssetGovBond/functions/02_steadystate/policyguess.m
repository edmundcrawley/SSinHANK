function [c_guess, inc]=policyguess(meshes,WW,RBRB,par,grid,N,W_fc,P_H,Profits_fc,Output)


jd_aux=P_H^1000;
jd_aux=jd_aux(1,:);

inc.labor   = par.tau.*WW.*meshes.h; % after tax labor income (labor disutility adjusted)
inc.money   = RBRB.*meshes.m;
inc.profits = sum((1-par.tau)*(N/par.H).*W_fc.*grid.h(1:end-1).*jd_aux(1:end-1))...
                  + (1-par.tau)*Profits_fc*par.profitshare*jd_aux(end)-(RBRB-1.0)*par.BtoY*Output; % lump sum transfer

%% Initial policy guesses:  Autarky policies as guess

% Consumption guess
c_guess = inc.labor + max(inc.money,0) + inc.profits;


end
