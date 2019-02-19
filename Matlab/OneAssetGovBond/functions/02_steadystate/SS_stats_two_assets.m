% Steady state statistics    
    %%    
    
    clear targets
    targets.ShareBorrower=sum((grid.m<0).*sum(joint_distr,2)');
    targets.B=grid.m*(sum(joint_distr,2));
    
    targets.BY=targets.B/Output;
    targets.Y=Output;
    BCaux_M=sum(sum(joint_distr,2),3);
    targets.m_bc=BCaux_M(1,:);
    targets.m_0=BCaux_M(grid.m==0);
   
    labortax =(1-par.tau).*W_fc.*N +(1-par.tau).*Profits_fc;
    par.gamma1=targets.B*(1-par.RB)+labortax(1);
    par.W=W_fc(1);
    par.PROFITS=Profits_fc(1);
    par.N=N;
    targets.GtoY=par.gamma1/Output;
    targets.T=labortax;
    
    %%   MPCs
[meshes.m,meshes.k,meshes.h] = ndgrid(grid.m,grid.k,grid.h);

NW=par.gamma/(1+par.gamma).*(par.N/par.H).*par.W;
WW=NW*ones(mpar.nm,mpar.nk,mpar.nh); %Wages
WW(:,:,end)=par.PROFITS*par.profitshare;
% MPC
WW_h=squeeze(WW(1,1,:));
WW_h_mesh=squeeze(WW(:,:,:).*meshes.h);

grid_h_aux=grid.h;

MPC_a_m = zeros(mpar.nm,mpar.nk,mpar.nh);
MPC_n_m = zeros(mpar.nm,mpar.nk,mpar.nh);
for kk=1:mpar.nk
    for hh=1:mpar.nh
        MPC_a_m(:,kk,hh)=gradient(squeeze(c_a_guess(:,kk,hh)))./gradient(grid.m)';
        MPC_n_m(:,kk,hh)=gradient(squeeze(c_n_guess(:,kk,hh)))./gradient(grid.m)';
    end
end

MPC_a_m=MPC_a_m.*(WW_h_mesh./c_a_guess);
MPC_n_m=MPC_n_m.*(WW_h_mesh./c_n_guess);

MPC_a_h = zeros(mpar.nm,mpar.nk,mpar.nh);
MPC_n_h = zeros(mpar.nm,mpar.nk,mpar.nh);
for mm=1:mpar.nm
    for kk=1:mpar.nk
        MPC_a_h(mm,kk,:)=gradient(squeeze(log(c_a_guess(mm,kk,:))))./gradient(log(WW_h'.*grid_h_aux))';
        MPC_n_h(mm,kk,:)=gradient(squeeze(log(c_n_guess(mm,kk,:))))./gradient(log(WW_h'.*grid_h_aux))';
    end
end

EMPC_h=joint_distr(:)'*(par.nu.*MPC_a_h(:)+(1-par.nu).*MPC_n_h(:));
EMPC_m=joint_distr(:)'*(par.nu.*MPC_a_m(:)+(1-par.nu).*MPC_n_m(:));

EMPC_a_h=joint_distr(:)'*MPC_a_h(:);
EMPC_a_m=joint_distr(:)'*MPC_a_m(:);

EMPC_n_h=joint_distr(:)'*MPC_n_h(:);
EMPC_n_m=joint_distr(:)'*MPC_n_m(:);

targets.Insurance_coeff=[1-EMPC_h 1-EMPC_m;
    1-EMPC_a_h 1-EMPC_a_m;
    1-EMPC_n_h 1-EMPC_n_m];