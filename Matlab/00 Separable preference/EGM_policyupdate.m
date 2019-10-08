function [c_star,m_star]=EGM_policyupdate(EVm,PIminus,RBminus,inc,meshes,grid,par,mpar,Hminus,Wminus,Profitminus)

%% EGM Step 1:
EMU = par.beta*reshape(EVm,[mpar.nm,mpar.nh]);
c_aux = 1./(EMU.^(1/par.xi));

mmin = min(grid.m);

% Calculate assets consistent with choices being (m')
% Calculate initial money position from the budget constraint,
% that leads to the optimal consumption choice
m_n_aux = (c_aux + meshes.m - c_aux.^(-par.xi/par.gamma).*(meshes.h/Hminus.*Wminus).^((1+par.gamma)/par.gamma))/(RBminus/PIminus);
m_n_aux(:,end) = (c_aux(:,end) + meshes.m(:,end) - inc.labor(:,end))/(RBminus/PIminus);
% m_n_aux = m_n_aux./(RBminus/PIminus+(m_n_aux<0)*par.borrwedge/PIminus);

% Identify binding constraints
binding_constraints = meshes.m < repmat(m_n_aux(1,:),[mpar.nm 1]);

[row,col]=find(binding_constraints);

options = optimset('Display','off');

cpbind  = zeros(mpar.nm,mpar.nh);

% cpbind = ((RBminus/PIminus)*meshes.m+mmin)/2+sqrt(((RBminus/PIminus)*meshes.m+mmin).^(2)+4*(meshes.h/Hminus*Wminus).^(1+par.varphi))/2;
% cpbind(:,end) = ((RBminus/PIminus)*meshes.m(:,end)+mmin+Profitminus*par.profitshare);

for i=1:length(row)
cpbind(row(i),col(i)) = fsolve( @(c) (RBminus/PIminus)*meshes.m(row(i),col(i))+...
           (c.^(-par.xi).*meshes.h(row(i),col(i))/Hminus.*Wminus).^(1/par.gamma).*meshes.h(row(i),col(i))/Hminus.*Wminus - c - mmin,c_aux(row(i),col(i)),options );
end

cpbind(:,end) = ((RBminus/PIminus)*meshes.m(:,end) - mmin+Profitminus*par.profitshare);



% Consumption when drawing assets m' to zero: Eat all Resources
% Resource = inc.labor + inc.money;

m_n_aux = reshape(m_n_aux,[mpar.nm mpar.nh]);
c_n_aux  = reshape(c_aux,[mpar.nm mpar.nh]);

% Interpolate grid.m and c_n_aux defined on m_n_aux over grid.m
% Check monotonicity of m_n_aux
if max(sum(abs(diff(sign(diff(m_n_aux))))))~=0
    warning('non monotone future liquid asset choice encountered')
end

c_star=zeros(mpar.nm,mpar.nh);
m_star=zeros(mpar.nm,mpar.nh);

for hh=1:mpar.nh
    Savings=griddedInterpolant(m_n_aux(:,hh),grid.m); % generate savings function a(s,a*)=a'
    m_star(:,hh)=Savings(grid.m); % Obtain m'(m,h) by Interpolation
    Consumption=griddedInterpolant(m_n_aux(:,hh),c_n_aux(:,hh)); % generate consumption function c(s,a*(s,a'))
    c_star(:,hh)=Consumption(grid.m);  % Obtain c(m,h) by interpolation (notice this is out of grid, used linear interpolation)
end

c_star(binding_constraints) = cpbind(binding_constraints);
% c_star(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_star(binding_constraints) = min(grid.m);

m_star(m_star>grid.m(end)) = grid.m(end);

end
