function [excess,c_new,L_new,Y,joint_distr,m_star,inc] = ...
           stationary_eqm(Guess,L_guess,jd_guess,crit,meshes,par,L,mpar,grid,P_H)

par.RB=Guess;

mc   =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;

% get dimensions of the grid
[nm,nh] = size(meshes.m);
mmin = meshes.m(1,1);
mmax = meshes.m(end,end);

distL    = crit+1;
iter_L = 0;

Y=sum(sum(L_guess.*jd_guess));
pf_guess=(1-mc)*Y;
% distance between two successive functions
% start with a distance above the convergence criterion
distC = crit+1;

iter = 0;

disp('Solving household problem by EGM')
tic

while (distL>crit)

% w0 = mc*par.alpha*(mc*(1-par.alpha)/(par.RB-1+par.delta))^((1-par.alpha)/par.alpha); % when capital is considered
w0 = mc; % labor share = 1 (labor is only one factor)
WW = w0.*L_guess;
WW(:,end) = pf_guess*par.profitshare;

% initial guess on future consumption 
% (consume asset income plus labor income from working l=1)
inc.labor = meshes.h/par.H.*WW;
inc.money = (par.RB)*meshes.m; % initial interest income
c_guess = inc.money+inc.labor; % guess on C_{t+1}

%%% iterate on the consumption decision rule
money_expense  = repmat(grid.m',[1 mpar.nh]);

distC = crit+1;

while (distC>crit)

 mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment

 mutil_c=par.RB.*mutil_c; %take return on money into account
 aux=reshape(permute(mutil_c,[2 1]),[mpar.nh mpar.nm]);
 % form expectations
 EMU_aux = par.beta*permute(reshape(P_H*aux,[mpar.nh mpar.nm]),[2 1]);

 c_aux = 1./(EMU_aux.^(1/par.xi)); % c0
     
 [c_new,m_star]=EGM_Step1_a(grid,inc,money_expense,c_aux,mpar,par,meshes,w0,mmin,pf_guess);
     
 m_star(m_star>grid.m(end)) = max(grid.m);

 % distance measure
 distC = max((abs(c_guess(:)-c_new(:))));
 % increase counter
 iter = iter+1;

 % update the guess on the consumption function
 c_guess = c_new;

end

% L  = @(c,h,w) invvp(up(c).*h.*w);

 L_new = L(c_new,[meshes.h(:,1:end-1)/par.H,zeros(mpar.nm,1)],w0); % Es do not supply labor

 distL = max((abs(L_guess(:)-L_new(:))));
 disp(distL);
 
 iter_L = iter_L+1;
 
 L_guess = L_new;

  [joint_distr]=JDiteration(m_star,P_H,mpar,grid);
 
 joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
 
 N = L_new(:)'*joint_distr(:);
%  K = grid.m*sum(joint_distr,2);
%  Y=N^(par.alpha)*K^(1-par.alpha);
 Y = N;
 pf_guess = (1-mc)*Y;
%  r_ss = mc*(1-par.alpha)*K^(-par.alpha)*N^(par.alpha) - par.delta;
%  par.RB=1+r_ss;

 
end

disp(([distC,distL]));
toc


fprintf('two loops for C and L, done. \n');

% market clearing condition
  AggregateSaving = m_star(:)'*joint_distr(:);
  grid.B=par.BtoY*Y; 
  excess(1)=(grid.B - AggregateSaving);


end
% 
