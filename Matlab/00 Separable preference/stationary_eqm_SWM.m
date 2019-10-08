function [cp0,c_update,L_new,joint_distr] = stationary_eqm_SWM(r0,crit,meshes,alpha,delta,rho,varphi,A0,C0,L,mpar,grid,P_H)

meshes_m = meshes.m;meshes_h = meshes.h;

% get dimensions of the grid
[nm,nh] = size(meshes_m);
mmin = meshes_m(1,1);
mmax = meshes_m(end,end);
% compute the wage according to marginal pricing and a Cobb-Douglas 
% production function
w0 = (1-alpha)*(alpha/(r0+delta))^(alpha/(1-alpha));

% initial guess on future consumption (consume asset income plus labor income from working h=1.
inc.labor = meshes_h*w0;
% inc.money = (1+r0)*meshes_m;
inc.money = (r0)*meshes_m;
cp0 = inc.money+inc.labor;

%%% iterate on the consumption decision rule

% distance between two successive functions
% start with a distance above the convergence criterion
dist_C    = crit+1;

% counter
iter    = 0;

dist_L    = crit+1;
iter_L = 0;
L_guess = ones(mpar.nm,mpar.nh);

fprintf('Inner loop, running... \n');

while (dist_L>crit)
% initial guess on future consumption (consume asset income plus labor income from working h=1.
inc.labor = meshes_h*w0.*L_guess;
inc.money = (1+r0)*meshes_m;
% inc.money = (r0)*meshes_m;
cp0 = inc.money+inc.labor;

%%% iterate on the consumption decision rule

% distance between two successive functions
% start with a distance above the convergence criterion
dist_C    = crit+1;

while (dist_C>crit)

% % current consumption level, cp0(anext,ynext) is the guess
% C0 = @(cp0,r) invup(beta*(1+r)*up(cp0)*pi');    
    
% derive current consumption, c0 == c_aux, calculated c_t from EE
c0 = C0(cp0,r0); % cp0 is the initial guess on future consumption

% % current asset level, c0 = C0(cp0(anext,ynext))
% A0 = @(anext,h,c0,r,w) 1/(1+r)*(c0+anext-L(c0,h,w).*h.*w);

% derive current assets, a0 is equivalent to m_star (initial money position)
a0 = A0(meshes_m,meshes_h,c0,r0,w0);

binding_constraints = meshes_m < repmat(a0(1,:),[mpar.nm 1 ]); % a0(1,:): the smallest current asset, thus, we think it's binding-constraint if saving is smaller than the smallest current asset

%%% update the guess for the consumption decision

% consumption decision rule for a binding borrowing constraint
% can be solved as a quadratic equation
% c^(2)-((1+r)a+b)c-(yw)^(1+2/3) = 0
% general formula and notation: ax^(2)+bx+c = 0
% x = (-b+sqrt(b^2-4ac))/2a
cpbind = ((1+r0)*meshes_m+mmin)/2+sqrt(((1+r0)*meshes_m+mmin).^(2)+4*(meshes_h*w0).^(1+varphi))/2;
% % slow alternative: rootfinding gives the same result
% cpbind  = zeros(nm,nh);
% options = optimset('Display','off');
% 
% for i=1:nh
% cpbind(:,i) = fsolve(@(c) (1+r0)*meshes_m(:,i)+L(c,grid_h(i),w0).*grid_h(i).*w0-c+mmin,cp0(:,i),options);
% end


% consumption for nonbinding borrowing constraint
cpnon = zeros(nm,nh);

% interpolation conditional on productivity realization instead of extrapolation use the highest value in a0
for i=1:nh
cpnon(:,i)  = interp1(a0(:,i),c0(:,i),meshes_m(:,i),'spline'); % meshes_m= anext, thus, cpnon=cnext (non-binding)
mpnon(:,i)  = interp1(a0(:,i),meshes_m(:,i),meshes_m(:,i),'spline');
end

mnext = mpnon;
beyond_mmax = mpnon > mmax;
mnext(beyond_mmax) = mmax;
mnext(binding_constraints) = mmin;


% % merge the two, cpbind and cpnon

% for i=1:nh % This holds only when the method is 'spline' in interp1. 'linear' gives NaN for the cpnext(binding_constraints)
% cpnext(:,i) = (meshes_m(:,i)>a0(1,i)).*cpnon(:,i)+(meshes_m(:,i)<=a0(1,i)).*cpbind(:,i);
% end

cpnext = cpnon;
cpnext(binding_constraints) = cpbind(binding_constraints);

m_star = reshape(a0,[mpar.nm mpar.nh]);
c_aux= reshape(c0,[mpar.nm mpar.nh]);

c_update = zeros(mpar.nm,mpar.nh); % c_update == cpnon
m_update = zeros(mpar.nm,mpar.nh);

for s=1:nh % produce consumption, saving policy function for each income grid (grid.s)
    M_choice=griddedInterpolant(m_star(:,s),grid.m,'spline'); % generate M choice function m(m*,k,h)=m' 
    m_update(:,s)=M_choice(grid.m); % Obtain m'(h,k,m) by Interpolation (notice this is out of grid, used linear interpolation)
    % This m_update is the interpolated (true) policy function : M(h,m)
    C_choice=griddedInterpolant(m_star(:,s),c_aux(:,s),'spline'); % generate consumption function c(s,a*(s,a'))
    c_update(:,s)=C_choice(grid.m);  % Obtain c(s,a) by interpolation (notice this is out of grid, used linear interpolation)
    % This C is the interpolated (true) policy function : C(s,a)           
end

c_update(binding_constraints) = cpbind(binding_constraints);
m_update(binding_constraints) = min(grid.m);
m_update(beyond_mmax) = max(grid.m);

% distance measure
dist_C = norm((cpnext-cp0)./cp0);

% display every 100th iteration
if mod(iter,100) == 1
    fprintf('Inner loop, iteration: %3i, Norm: %2.6f \n',[iter,dist_C]);
end

% increase counter
iter = iter+1;

% update the guess on the consumption function
cp0 = cpnext;

end

 L_new = L(c_aux,meshes.h,w0); % Es do not supply labor

 dist_L = max((abs(L_guess(:)-L_new(:))));
 disp(dist_L);
 
 iter_L = iter_L+1;
 
 L_guess = L_new;


end
fprintf('Inner loop, done. \n');


  [joint_distr]=JDiteration(m_update,P_H,mpar,grid);

  joint_distr=reshape(joint_distr,[mpar.nm mpar.nh]);
  

end

