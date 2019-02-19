function [c_guess, inc]=policyguess(meshes,WW,RBRB,par,mpar)

inc.labor   = par.tau.*WW.*meshes.h;
inc.money   = RBRB.*meshes.m;
inc.profits = 0;%lump sum profits

%% Initial policy guesses:  Autarky policies as guess

% Consumption guess
c_guess = inc.labor + max(inc.money,0) + inc.profits;


end
