function [excess,c_guess,c_update,c_aux,K,N,r_ss,Profits_fc,joint_distr,w0] = stationary_eqm_guess1(Guess,crit,meshes,par,A0,C0,L,mpar,grid,P_H)

par_RB = Guess(1);
Profits_fc = Guess(2);

mc   =  par.mu - (par.beta * log(par.PI) - log(par.PI))/par.kappa;

% get dimensions of the grid
[nm,nh] = size(meshes.m);
mmin = meshes.m(1,1);
mmax = meshes.m(end,end);

distL    = crit+1;
iter_L = 1;
L_guess = ones(mpar.nm,mpar.nh);
L_guess(:,end)=0;

% distance between two successive functions
% start with a distance above the convergence criterion
distC = crit+1;

iter = 1;

disp('Solving household problem by EGM')
tic

while (distL>crit)

w0 = mc*par.alpha*(mc*(1-par.alpha)/(par_RB-1+par.delta))^((1-par.alpha)/par.alpha);
WW = w0.*L_guess;
WW(:,end) = Profits_fc*par.profitshare;

% initial guess on future consumption 
% (consume asset income plus labor income from working l=1)
inc_labor = meshes.h.*WW;
inc_money = (par_RB-1)*meshes.m; % initial interest income
c_guess = inc_money+inc_labor; % guess on C_{t+1}

%%% iterate on the consumption decision rule



while (distC>crit)

     mutil_c = 1./(c_guess.^par.xi); % marginal utility at consumption policy no adjustment

     mutil_c=par_RB.*mutil_c; %take return on money into account
     aux=reshape(permute(mutil_c,[2 1]),[mpar.nh mpar.nm]);
     % form expectations
     EMU_aux = par.beta*permute(reshape(P_H*aux,[mpar.nh mpar.nm]),[2 1]);

     c_aux = 1./(EMU_aux.^(1/par.xi)); % c0
     m_star = (c_aux + meshes.m - inc_labor)./par.RB; % a0

     binding_constraints = meshes.m < repmat(m_star(1,:),[mpar.nm 1 ]); % a0(1,:): the smallest current asset, thus, we think it's binding-constraint if saving is smaller than the smallest current asset

     Resource = inc_labor  + inc_money;

     %%% update the guess for the consumption decision

     % consumption decision rule for a binding borrowing constraint
     % can be solved as a quadratic equation
     % c^(2)-((1+r)a+b)c-(yw)^(1+2/3) = 0
     % general formula and notation: ax^(2)+bx+c = 0
     % x = (-b+sqrt(b^2-4ac))/2a
%      cpbind = ((par_RB)*meshes.m+mmin)/2+sqrt(((par_RB)*meshes.m+mmin).^(2)+4*(meshes.h*w0).^(1+par.varphi))/2;
%      cpbind(:,end) = par_RB*meshes.m(:,end) + mmin + Profits_fc*par.profitshare;
     % % slow alternative: rootfinding gives the same result
     % cpbind  = zeros(nm,nh);
     % options = optimset('Display','off');
     % 
     % for i=1:nh
     % cpbind(:,i) = fsolve(@(c) (1+r0)*meshes.m(:,i)+L(c,grid.h(i),w0).*grid.h(i).*w0-c+mmin,cp0(:,i),options);
     % end


     % consumption for nonbinding borrowing constraint
     cpnon = zeros(nm,nh);

     % interpolation conditional on productivity realization instead of extrapolation use the highest value in a0
     for i=1:nh
     cpnon(:,i)  = interp1(m_star(:,i),c_aux(:,i),meshes.m(:,i),'spline'); % meshes.m= anext, thus, cpnon=cnext (non-binding)
     mpnon(:,i)  = interp1(m_star(:,i),meshes.m(:,i),meshes.m(:,i),'spline');
     end

     mnext = mpnon;
     beyond_mmax = mpnon > mmax;
     mnext(beyond_mmax) = mmax;
     mnext(binding_constraints) = mmin;


     % % merge the two, cpbind and cpnon

     % for i=1:nh % This holds only when the method is 'spline' in interp1. 'linear' gives NaN for the cpnext(binding_constraints)
     % cpnext(:,i) = (meshes.m(:,i)>a0(1,i)).*cpnon(:,i)+(meshes.m(:,i)<=a0(1,i)).*cpbind(:,i);     
     % end

     cpnext = cpnon;
     cpnext(binding_constraints) = Resource(binding_constraints)-grid.m(1);
%      cpnext(binding_constraints) = cpbind(binding_constraints);

     m_star = reshape(m_star,[mpar.nm mpar.nh]);
     % c_aux= reshape(c0,[mpar.nm mpar.nh]);

     c_update = zeros(mpar.nm,mpar.nh); % c_update == cpnon
     m_update = zeros(mpar.nm,mpar.nh);

     for s=1:nh % produce consumption, saving policy function for each income grid (grid.s)
         M_choice=griddedInterpolant(m_star(:,s),grid.m,'spline'); 
         m_update(:,s)=M_choice(grid.m); 
         C_choice=griddedInterpolant(m_star(:,s),c_aux(:,s),'spline'); 
         c_update(:,s)=C_choice(grid.m);  
         
     end

     c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
%      c_update(binding_constraints) = cpbind(binding_constraints);

     m_update(binding_constraints) = min(grid.m);
     m_update(beyond_mmax) = max(grid.m);

     % distance measure
%      distC = max((abs(cpnext(:)-c_guess(:))));
     distC = max((abs(c_update(:)-c_guess(:))));
     % increase counter
     iter = iter+1;

     % update the guess on the consumption function
     c_guess = cpnext;
%      disp(([distC]));
end

% L  = @(c,h,w) invvp(up(c).*h.*w);

L_new = L(c_aux,[meshes.h(:,1:end-1),zeros(mpar.nm,1)],w0); % Es do not supply labor

distL = max((abs(L_new(:)-L_guess(:))));
disp(distL);
iter_L = iter_L+1;
L_guess = L_new;

end

disp(([distC,distL]));
toc

fprintf('Outer loop, done. \n');

  LS_ind = L(c_aux,grid.h,w0); %LS_ind = L(c_aux,grid.h,w);
  LS_ind(:,end)=0;
  [joint_distr]=JDiteration(m_star,P_H,mpar,grid);

  joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
  N = LS_ind(:)'*joint_distr(:);
  K = grid.m*sum(joint_distr,2);
  Y=N^(par.alpha)*K^(1-par.alpha);
  Profits_fc_ep = 1*(1-mc)*Y - Y.*(1/(1-par.mu))./par.kappa./2 .*log(par.PI).^2;
  r_ss = mc*(1-par.alpha)*K^(-par.alpha)*N^(par.alpha) - par.delta;
  
  excess(1) = par_RB - (1+r_ss);
  excess(2) = Profits_fc - Profits_fc_ep;

end

