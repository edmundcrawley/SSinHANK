function [c_update,m_update]=EGM_Step1_a(grid,inc,money_expense,c_aux,mpar,par,meshes,w0,mmin,pf_guess)

cpbind = ((par.RB)*meshes.m+mmin)/2+sqrt(((par.RB)*meshes.m+mmin).^(2)+4*(meshes.h/par.H*w0).^(1+par.varphi))/2;

% cpbind  = zeros(mpar.nm,mpar.nh);
% options = optimset('Display','off');
% 
% for i=1:mpar.nh-1
% cpbind(:,i) = fsolve(@(c) (par.RB)*meshes.m(:,i)+...
%                (c.^(-par.xi).*meshes.h(:,i)/par.H.*w0).^(1/par.gamma).*meshes.h(:,i)/par.H.*w0-c+mmin,c_aux(:,i),options);
% end

cpbind(:,end) = ((par.RB)*meshes.m(:,end)+mmin+pf_guess*par.profitshare);

m_star = (c_aux + money_expense - inc.labor)./par.RB; % a0

binding_constraints = money_expense < repmat(m_star(1,:),[mpar.nm 1 ]); % a0(1,:): the smallest current asset, thus, we think it's binding-constraint if saving is smaller than the smallest current asset

% c_aux(binding_constraints)=cpbind(binding_constraints);
% 
% m_star = (c_aux + money_expense - inc.labor)./par.RB; % a0

% Resource = inc.labor  + inc.money;

m_star = reshape(m_star,[mpar.nm mpar.nh]);
     
c_update = zeros(mpar.nm,mpar.nh); % c_update == cpnon
m_update = zeros(mpar.nm,mpar.nh);

for s=1:mpar.nh % produce consumption, saving policy function for each income grid (grid.s)
%    M_choice=griddedInterpolant(m_star(:,s),grid.m,'spline'); 
   M_choice=griddedInterpolant(m_star(:,s),grid.m); 
   m_update(:,s)=M_choice(grid.m); 
   C_choice=griddedInterpolant(m_star(:,s),c_aux(:,s)); 
   c_update(:,s)=C_choice(grid.m);  
         
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new handling of no borrowing
% binding_constraints = c_update>Resource;
% c_update(binding_constraints) = Resource(binding_constraints);
% m_update(binding_constraints) = min(grid.m);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


c_update(binding_constraints) = cpbind(binding_constraints);
% c_update(binding_constraints) = Resource(binding_constraints)-grid.m(1);
m_update(binding_constraints) = min(grid.m);

end
